/*----------------------------------------------------------------------
 *                           Program Bdd-bin-dse
 *----------------------------------------------------------------------
 *
 *  In this test program, one processor reads the whole matrix
 *  from file and broadcast to other processors. The matrix is
 *  assumed to be in binary format. Matrix graph is then  partitioned
 *  using  DSE, a simple partitioning routine, and scatters the local
 *  matrices to each processor. Once these submatrices are received each
 *  processor solves the problem using preconditioned FGMRES
 *  preconditioned with :
 *                       BJ, RAS, SCHUR
 *--------------------------------------------------------------------*/

#if defined(__ICC)
#include <mathimf.h>
#else
#include <math.h>
#endif
#include "aux.h"

int main(int argc, char *argv[])
{

    /* declarations related to binary format for reading the bin
     matri. Second part is related to I/O parameters */
    char mname[MAX_MAT][MAX_LINE], key[MAX_LINE];
    int nrhs, mat;
    long int m;    // = header(1); //row;
    long int n;    // = header(2); //column;
    long int nnz;    // = header(3);%nnz
    long int header; // header to store file type, mat or rhs;
    long int numread;

    int myid, ierr, i, nloc;
    /* memory usage of the preconditioning matrix */
    double ratio, Bdensity;
    /* working array for reading matrix */
    double norm, res1, tpc, ttol;
    FLOAT *a, *rhstmp;
    int *ja, *nnzptr, *gia;
    int npro,its, *im, *bim;//bim is for block case
    fprm prm;
    FILE *fp=NULL, *fout=NULL, *binfile=NULL;
    char *name, *iluname, *curmat, *currhs, *curname = NULL, buf[40];
    int *nzding = NULL;
    int nBlock, *nB = NULL, *perm = NULL;
    double tib1, tib2, tib3, tib4 = 0.0, blocksize = 0.0;

    /*-------------------- variables related to dse partitioning */
    int *riord, *dom = NULL, *idom = NULL, *mask, *jwk, *link;

    /*-------------------- variables related to pARMS */
    parms_Map       map;
    FLOAT       *x, *y, *rhs, *resvec;
    parms_Mat       A;
    parms_PC        pc;
    parms_Solver    solver;
    parms_Timer     tm;


    /*-------------------- initialize MPI environment */
    MPI_Init(&argc, &argv);//initilization
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);//return the pid
    MPI_Comm_size(MPI_COMM_WORLD, &npro);//return the number of the processors

    nrhs = 0;

    /*-------------------- read matrix name from input file */
    prm = malloc(sizeof(*prm));
    if (prm == NULL) {
        fprintf(stderr, "cannot allocate memory for prm\n");
        MPI_Abort(MPI_COMM_WORLD, 66);
    }

    if (argc >= 2) {
        read_param(argv[1], mname, prm);
    }
    else if (argc == 1){
        read_param("inputs", mname, prm);
    }
    parms_TimerCreate(&tm);

    /* --- Begin loop over matrices ----*/
    mat = 0;
    while(mname[mat][1] != '#'){
        a = NULL; ja = NULL; nnzptr = NULL; rhstmp = NULL;
        curmat = strtok(mname[mat], " ");
        if (curmat == NULL)
            continue;

        currhs = strtok(NULL, " ");

        if (myid == 0)
            printf("Reading matrix %s\n", curmat);

        if ((binfile = fopen(curmat, "r")) == NULL) {
            fprintf(stderr, "Error opening matrix file\n");
            fprintf(stderr, "filename = %s\n", curmat);
            MPI_Finalize();
            exit(1);
        }

        numread = fread(&header, sizeof(int), 1, binfile);
        header = bswap_32(header);
        printf("header is %ld \n", header );

        if(header == 0) {
            fprintf(stderr, "File does not have that many items");
            MPI_Finalize();
            exit(1);
        }

        numread = fread(&m, sizeof(int), 1, binfile);
        m = bswap_32(m);//row
        printf("m is %ld \n", m );

        numread = fread(&n, sizeof(int), 1, binfile);
        n = bswap_32(n);//column
        printf("n is %ld \n", n );

        if(n != m) {
            fprintf(stderr, "Matrix is not square");
            MPI_Finalize();
            exit(1);
        }

        numread = fread(&nnz, sizeof(int), 1, binfile);
        nnz = bswap_32(nnz);// nnz

        nnzptr = malloc(m*sizeof(*nnzptr)); //nnz of each row
        numread = fread(nnzptr, sizeof(int), m, binfile);

        if (numread != m){
            fprintf(stderr," error in fread");
            MPI_Finalize();
            exit(1);
        }

        gia = malloc((n+1)*sizeof(*gia));
        gia[0] = 0;

        long int sum_nz = 0;
        for (i = 0; i < m; ++i){
            nnzptr[i] = bswap_32(nnzptr[i]);
            sum_nz += nnzptr[i];
            gia[i+1] = sum_nz;
        }

        if (myid == 0){
            if(sum_nz != nnz){
                fprintf(stderr," No-Nonzeros sum-rowlengths do not match %ld %ld", nnz, sum_nz);
                MPI_Finalize();
                exit(1);
            }

            ja = malloc(nnz*sizeof(*ja)); //column indeces of all non-zero entries
            numread = fread(ja, sizeof(int), nnz, binfile);
            if (numread != nnz){
                fprintf(stderr," error in fread");
                MPI_Finalize();
                exit(1);
            }

            for (i = 0; i < nnz; ++i)
                ja[i] = bswap_32(ja[i]);//swap between big endian and small endian

            a = malloc(nnz*sizeof(*a)); //values of all non-zero entries
            if (a == NULL) {
                fprintf(stderr, "cannot allocate memory for a\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }

            numread = fread(a, sizeof(double), nnz, binfile);
            if (numread != nnz){
                fprintf(stderr," error in fread");
                return 1;
            }

            for (i = 0; i < nnz; ++i)
                a[i] = byteswap_double(a[i]);

            fseek(binfile, 0, 0);

            if(myid == 0){
                if(argc == 3) {
                    if (NULL == (fp = fopen(argv[2], "w"))) {
                        fprintf(stderr, "cannot open file %s\n", argv[2]);
                        exit(1);
                    }
                }
                else {
                    fp = stdout;
                }
                curname = strtok(curmat, "/");
                while (curname != NULL)
                {
                    strcpy(key, curname);
                    curname = strtok(NULL, "/");
                }
                curname = strtok(key, ".");
                fprintf(fp, "\nMatrix %d: %s \n",(mat+1), curname);
                fprintf(fp, "n = %ld, nnz = %ld\n", n, nnz);
            }
            parms_TimerReset(tm);

            csptr csmat = NULL;

            csmat = malloc(sizeof(*csmat));

            bincols2csptr(n, nnzptr, ja, a, csmat);

            free(a);
            free(ja);

            vbsptr vbmat = NULL;

            printf("prm->eps = %f\n",prm->eps);

            if (prm->cosine)
                ierr = init_blocks( csmat, &nBlock, &nB, &perm, prm->eps);
            else
                ierr = init_blocks_density( csmat, &nBlock, &nB, &perm, prm->eps);
            printf("prm->cosine = %d\n",prm->cosine);

            tib1 =  parms_TimerGet(tm);
            printf("\ntime on init=%f\n",tib1);
            if(ierr != 0) {
                fprintf(stderr, "*** in init_blocks ierr != 0 ***\n");
                MPI_Finalize();
                exit(1);
            }

            if( dpermC( csmat, perm ) != 0 ) {
                fprintf( stderr, "*** dpermC error ***\n" );
                MPI_Finalize();
                exit(1);
            }
            tib2 =  parms_TimerGet(tm);
            printf("\ntime on dpermC=%f\n",tib2-tib1);
            /*-------------------- convert to block matrix. */
            vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
            ierr = csrvbsrC_new( 1, nBlock, nB, csmat, vbmat );
            tib3 =  parms_TimerGet(tm);
            printf("\ntime on csrvbsrC_new=%f\n",tib3-tib2);


            blocksize = (double)csmat->n / (double)nBlock;
            Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;
            if (myid == 0)
                printf("\n Bsize=%-7f, Bdensity=%-7f\n",blocksize, Bdensity);




            int nbb, *bia, *bja;
            BData *ba;
            bia = (int*)malloc((vbmat->n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
            bja = (int*)malloc(nnz*sizeof(int));
            ba = (BData*)malloc(nnz*sizeof(BData));
            vbsptr2colunms(vbmat, &nbb, bia, bja, ba);//int vbsptr2colunms(vbsptr mat, int *n, int *ia, int *ja)

            bja = (int*)realloc(bja, nbb*sizeof(int));

            tib4 =  parms_TimerGet(tm);
            printf("\ntime on csrvbsrC_new=%f, time on whole blocking process = %f\n",tib4-tib3, tib4);
            /*--------------------get the elapsed time spent on creating PC */

            idom = malloc((npro+1)*sizeof(*idom));
            if (idom == NULL) {
                fprintf(stderr, "cannot allocate memory for idom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            dom = malloc(nBlock*sizeof(*dom));
            if (dom == NULL) {
                fprintf(stderr, "cannot allocate memory for dom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }

            if (npro == 1) {
                for (i = 0; i < nBlock; i++) {
                    dom[i] = i+1;
                }
                idom[0] = 1;
                idom[1] = nBlock+1;
            }
            else {
                riord = malloc(nBlock*sizeof(*riord));
                if (riord == NULL) {
                    fprintf(stderr, "cannot allocate memory for riord\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                mask = malloc(nBlock*sizeof(*mask));
                if (mask == NULL) {
                    fprintf(stderr, "cannot allocate memory for mask\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                jwk = malloc(2*nBlock*sizeof(*jwk));
                if (jwk == NULL) {
                    fprintf(stderr, "cannot allocate memory for jwk\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                link = malloc(nBlock*sizeof(*link));
                if (link == NULL) {
                    fprintf(stderr, "cannot allocate memory for link\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                dse_(&nBlock, bja, bia, &npro, riord, dom, idom, mask, jwk, link);
                free(riord);
                free(mask);
                free(jwk);
                free(link);
            }
            free(ba);
            free(bja);
            free(bia);


            nzding = (int *) calloc(npro, sizeof(int));
            int i1, i3;
            int i2 = 0, i4;
            for (i1 = 0; i1 < nBlock; ++i1)
            {
                if (i1 == idom[i2+1]-1)
                    i2++;
                i4 = dom[i1]-1;
                for (i3 = vbmat->bsz[i4]; i3 < vbmat->bsz[i4+1]; ++i3)
                    nzding[i2] += csmat->nnzrow[i3];
            }

            cleanCS( csmat );
            cleanVBMat( vbmat );
        }


        MPI_Bcast(&nBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (myid != 0)
        {
            nB = (int *) malloc(nBlock * sizeof(int));
            if (nB == NULL) {
                fprintf(stderr, "cannot allocate memory for nB\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            perm = (int *) malloc(n * sizeof(int));
            if (perm == NULL) {
                fprintf(stderr, "cannot allocate memory for perm\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            idom = malloc((npro+1)*sizeof(*idom));
            if (idom == NULL) {
                fprintf(stderr, "cannot allocate memory for idom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            dom = malloc(nBlock*sizeof(*dom));
            if (dom == NULL) {
                fprintf(stderr, "cannot allocate memory for dom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            nzding = (int *) malloc(npro * sizeof(int));
        }
        MPI_Bcast(nB, nBlock, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(perm, n, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(idom, npro+1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(dom, nBlock, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(nzding, npro, MPI_INT, 0, MPI_COMM_WORLD);

        /*-------------------- Create map object */
        parms_MapCreateFromPtr(&map, n, dom, idom, MPI_COMM_WORLD, 1, NONINTERLACED);

        parms_Map_Assign_blockstructure(map, nB);


        nloc = parms_MapGetLocalSize(map);

        a = malloc(nzding[myid]*sizeof(*a));
        if (a == NULL) {
            fprintf(stderr, "cannot allocate memory for a\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        ja = malloc(nzding[myid]*sizeof(*ja));
        if (ja == NULL) {
            fprintf(stderr, "cannot allocate memory for ja\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        ierr = local_read_bin_data_b(binfile, m, n, nnz, nnzptr, ja, a, idom, dom, perm, nB, nBlock, gia);
        fclose(binfile);


        free(dom);
        free(idom);
        free(gia);

        if(ierr != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "cannot read matrix\n");
            fprintf(stderr, "filename = %s\n", mname[mat]);
            MPI_Finalize();
            exit(1);
        }


        rhstmp = malloc(n*sizeof(*rhstmp));
        if (rhstmp == NULL) {
            fprintf(stderr, "cannot allocate memory for rhstmp\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        if (currhs != NULL)
        {
            if ((binfile = fopen(currhs, "r")) == NULL) {
                fprintf(stderr, "Error opening rhs file\n");
                MPI_Finalize();
                exit(1);
            }

            numread = fread(&header, sizeof(int), 1, binfile);
            header = bswap_32(header);
            printf("header is %ld \n", header );

            numread = fread(&m, sizeof(int), 1, binfile);
            m = bswap_32(m);//row

            if(m != n) {
                fprintf(stderr, "rhs's dim is different from mat's");
                MPI_Finalize();
                exit(1);
            }
            nrhs = 1;
            numread = fread(rhstmp, sizeof(double), m, binfile);
            if (numread != m){
                fprintf(stderr," error in fread");
                MPI_Finalize();
                exit(1);
            }

            for (i = 0; i < m; ++i)
                rhstmp[i] = byteswap_double(rhstmp[i]);//need to be optimized

            fclose(binfile);
        }
        else nrhs = 0;

        csptr csmat = NULL;

        csmat = malloc(sizeof(*csmat));


        bincols2csptr(n, nnzptr, ja, a, csmat);


        free(nnzptr);
        free(a);
        free(ja);

        vbsptr vbmat = NULL;
        /*--------------------create a timer */

        if( dpermC( csmat, perm ) != 0 ) {
            fprintf( stderr, "*** dpermC error ***\n" );
            MPI_Finalize();
            exit(1);
        }
        /*-------------------- convert to block matrix. */
        vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
        ierr = csrvbsrC_new( 1, nBlock, nB, csmat, vbmat );


        FLOAT *rhstmpp;
        rhstmpp = (FLOAT*)malloc(n*sizeof(FLOAT));


        for( i = 0; i < n; i++ )
            rhstmpp[perm[i]] = rhstmp[i];

        int nbb, *bia, *bja;
        BData *ba;
        bia = (int*)malloc((vbmat->n+1)*sizeof(int));
        bja = (int*)malloc(nzding[myid]*sizeof(int));
        ba = (BData*)malloc(nzding[myid]*sizeof(BData));
        vbsptr2colunms(vbmat, &nbb, bia, bja, ba);
        bja = (int*)realloc(bja, nbb*sizeof(int));


        free(nzding);

        /*-------------------- create a distributed matrix based on map */
        parms_MatCreate(&A, map);
        /*-------------------- Insert values into A */
        /* Initialize global indices array */
        im = (int *)malloc(n*sizeof(int));
        for (i = 0; i < n; i++) {
            im[i] = i+1;
        }

        bim = (int *)malloc(nBlock*sizeof(int));
        for (i = 0; i < nBlock; i++) {
            bim[i] = i+1;
        }
        parms_MatSetValues_b(A, nBlock, bim, bia, bja, ba, INSERT);

        free(ba);
        free(bja);
        free(bia);

        int llsize;
        //printf("llsizetest=%d\n", llsize);
        llsize = parms_MapGetLocallSize(map);
        printf("llsize=%d, nloc = %d\n", llsize, nloc);
        /*-------------------- Create distributed vectors based on map */
        x = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        rhs = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        y = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        resvec = (FLOAT *)malloc(llsize*sizeof(FLOAT));


        /*-------------------- Setup the matrix and communication structure */
        parms_Mat_B_Setup(A);
        printf("\n npro=%d\n",npro);

        /*-------------------- Copy rhs or Setup artifical right-hand-side vector */

        if(!nrhs)
        {
            for(i=0; i<llsize; i++)
            {
                x[i] = 1.0;
            }
            parms_MatVec(A, x, rhs);
        }
        else
        {
            parms_VecSetValues_b(rhs, nBlock, bim, rhstmpp, INSERT, map);
        }

        free(rhstmp);
        free(rhstmpp);
        free(im);

        free(bim);

        /*--------------------Setup initial guess to a vector of all-0. */
        for(i=0; i<llsize; i++)
        {
            x[i] = 0.0;
        }


        /*--------------------Get 2-norm of initial residual */

        parms_MatVec(A, x, resvec);

        for(i=0; i<llsize; i++)
        {
            resvec[i] = resvec[i] - rhs[i];
        }

        parms_VecGetNorm2_b(resvec, &norm, map);

        free(resvec);
        /*--------------------Create preconditioner based on the matrix A. */
        parms_PCCreate_b(&pc, A);

        /*--------------------set parameters for pc */
        set_pc_params_b(pc, prm);


        /*--------------------reset the timer */
        parms_TimerReset(tm);
        parms_PCSetup_b(pc);


        /*--------------------get the elapsed time spent on creating PC */
        tpc =  parms_TimerGet(tm);

        printf("The time for pc setup %f in proc %d\n", tpc, myid);


        /*--------------------get the ratio of the number of nonzero entries in the
    preconditioning matrix to that in the original matrix */
        parms_PCGetRatio(pc, &ratio);//block version is not ready
        /*--------------------pause the timer */
        parms_TimerPause(tm);



        /*--------------------Create a solver based on A and pc */

        parms_SolverCreate(&solver, A, pc);
        /*--------------------Set the solver type */
        parms_SolverSetType_b(solver, SOLFGMRES);



        /*--------------------set parameters for solver */
        set_solver_params(solver, prm);



        /*--------------------set up solver -- no longer needed - DOK */

        /*--------------------restart the timer */
        parms_TimerRestart(tm);

        /*--------------------Solve the linear equation */
        parms_SolverApply_b(solver, rhs, x);


        /*--------------------get total time spent on creating the pc and solving the linear
     system */
        ttol = parms_TimerGet(tm);

        printf("The time for solving  %f in proc %d\n", ttol-tpc, myid);
        printf("The total time cost  %f in proc %d\n", ttol, myid);

        /*--------------------Get the number of iterations */
        its = parms_SolverGetIts(solver);

        /*--------------------Compute the residual error  */
        parms_MatVec(A, x, y);
        for(i=0; i<llsize; i++)
        {
            y[i] = rhs[i] - y[i];
        }
        parms_VecGetNorm2_b(y, &res1, map);
        /*--------------------processor 0 outputs the result */
        if (myid == 0) {
            parms_PCGetName(pc, &name);
            parms_PCILUGetName(pc, &iluname);
            fprintf(fp, "The blocking information Bsize=%-7f, Bdensity=%-7f\n", blocksize, Bdensity);
            fprintf(fp, "The preconditioner  %9s %s\n", " ", name);
            fprintf(fp, "The local preconditioner %4s %s\n", " ", iluname);
            fprintf(fp, "The memory usage %12s %-4.2f\n", "=", ratio);
            fprintf(fp, "The number of processors %4s %-4d\n", "=", npro);
            fprintf(fp, "The number of iterations %4s %-4d\n", "=", its);
            sprintf(buf, "%8.2fs", tpc);
            fprintf(fp, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol-tpc);
            fprintf(fp, "The solving time %12s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol);
            fprintf(fp, "The total time %14s %-s\n", "=", strtok(buf, " "));
            fprintf(fp, "The initial residual norm %3s %-8.2e\n", "=",norm);
            fprintf(fp, "The final residual norm %5s %-8.2e\n","=", res1);

            if(fp != stdout)
                fclose(fp);
            /* --- Write output to file ---- */
            fout = fopen("output.txt", "aw");
            fprintf(fout, "\nMatrix %d: %s \n",(mat+1), curname);
            fprintf(fout, "n = %ld, nnz = %ld\n", n, nnz);
            fprintf(fout, "The blocking information Bsize=%-7f, Bdensity=%-7f\n", blocksize, Bdensity);
            fprintf(fout, "The preconditioner  %9s %s\n", " ", name);
            fprintf(fout, "The local preconditioner %4s %s\n", " ", iluname);
            fprintf(fout, "The memory usage %12s %-4.2f\n", "=", ratio);
            fprintf(fout, "The number of processors %4s %-4d\n", "=", npro);
            fprintf(fout, "The number of iterations %4s %-4d\n", "=", its);
            sprintf(buf, "%8.2fs", tpc);
            fprintf(fout, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol-tpc);
            fprintf(fout, "The solving time %12s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol);
            fprintf(fout, "The total time %14s %-s\n", "=", strtok(buf, " "));
            fprintf(fout, "The initial residual norm %3s %-8.2e\n", "=",norm);
            fprintf(fout, "The final residual norm %5s %-8.2e\n","=", res1);
            fclose(fout);
        }

        fflush(stdout);
        MPI_Barrier(MPI_COMM_WORLD);

        /*--------------------Free memories */
        cleanCS( csmat );
        cleanVBMat( vbmat );


        free(x);
        free(y);
        free(rhs);
        free(perm);
        parms_MatFree_b(&A);
        parms_MapFree(&map);
        parms_PCFree_b(&pc);
        parms_SolverFree_b(&solver);
        parms_TimerFree(&tm);

        /*----Goto next matrix ---*/
        //next_mat://for block structure detection
        mat++;

    }

    /*--------------Free prm--------------*/
    free(prm);
    /*--------------------Exit MPI environment */
    MPI_Finalize();
    return 0;
}

