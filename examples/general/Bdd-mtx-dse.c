//#define _GNU_SOURCE
/*----------------------------------------------------------------------
 *                           Program dd-HB-dse
 *----------------------------------------------------------------------
 *
 *  In this test program, each processor reads the whole matrix
 *  from file. The matrix is assumed to be in Harwell-Boeing format.
 *  Matrix graph is then  partitioned  using  DSE, a simple partitioning
 *  routine, and scatters the local matrices to each processor. Once
 *  these submatrices are received each processor solves the problem
 *  using preconditioned FGMRES preconditioned with :
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

    /* declarations related to matrix market format for reading the mtx
     matrices. Second part is related to I/O parameters */
    char mname[MAX_MAT][MAX_LINE], key[MAX_LINE];//, type[3];//guesol[2], title[72],
    int nrhs, nc, n, nnz, mat;//tmp0,tmp2, tmp3, tmp,
    int myid, ierr, i, nloc;
    /* memory usage of the preconditioning matrix */
    double ratio, Bdensity = 0.0;
    /* working array for reading matrix */
    double norm, res1, tpc, ttol;
    FLOAT *a, *rhstmp;
    int *ja, *ia;
    int npro,its, *im, *bim;//bim is for block case
    fprm prm;
    FILE *fp=NULL, *fout=NULL, *mtxfile=NULL;
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
        read_param(argv[1],mname, prm);
    }
    else if (argc == 1){
        read_param("inputs",mname, prm);
    }


    /* --- Begin loop over matrices ----*/
    mat = 0;
    while(mname[mat][1] != '#'){
        a = NULL; ja = NULL; ia = NULL; rhstmp = NULL;
        curmat = strtok(mname[mat], " ");
        if (curmat == NULL)
            continue;
        currhs = strtok(NULL, " ");
        parms_TimerCreate(&tm);
        parms_TimerReset(tm);

        if (myid == 0)
            printf("Reading matrix %s\n", curmat);

        if ((mtxfile = fopen(curmat, "r")) == NULL) {
            fprintf(stderr, "Error opening matrix file\n");
            fprintf(stderr, "filename = %s\n", curmat);
            MPI_Finalize();
            exit(1);
        }

        MM_typecode matcode;
        ierr = mm_read_banner(mtxfile, &matcode);
        if(ierr != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "cannot read banner\n");
            fprintf(stderr, "filename = %s\n", curmat);
            MPI_Finalize();
            exit(1);
        }

        if ((ierr = mm_read_mtx_crd_size(mtxfile, &n, &nc, &nnz, matcode)) != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "Error reading matrix dimensions\n");
            MPI_Finalize();
            exit(1);
        }

        if (myid == 0)
        {
            a = malloc(nnz*sizeof(*a));
            if (a == NULL) {
                fprintf(stderr, "cannot allocate memory for a\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            ja = malloc(nnz*sizeof(*ja));
            if (ja == NULL) {
                fprintf(stderr, "cannot allocate memory for ja\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            ia = malloc(nnz*sizeof(*ia));
            if (ia == NULL) {
                fprintf(stderr, "cannot allocate memory for ia\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }


            ierr = mm_read_mtx_crd_data(mtxfile, n, nc, nnz, ia, ja, a, matcode);
            fclose(mtxfile);

            if ((mtxfile = fopen(curmat, "r")) == NULL) {
                fprintf(stderr, "Error opening matrix file\n");
                fprintf(stderr, "filename = %s\n", curmat);
                MPI_Finalize();
                exit(1);
            }

            MM_typecode matcode;
            ierr = mm_read_banner(mtxfile, &matcode);
            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read banner\n");
                fprintf(stderr, "filename = %s\n", curmat);
                MPI_Finalize();
                exit(1);
            }

            if ((ierr = mm_read_mtx_crd_size(mtxfile, &n, &nc, &nnz, matcode)) != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "Error reading matrix dimensions\n");
                MPI_Finalize();
                exit(1);
            }

            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read matrix\n");
                fprintf(stderr, "filename = %s\n", mname[mat]);
                MPI_Finalize();
                exit(1);
            }

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
                fprintf(fp, "n = %d, nnz = %d\n", n, nnz);
            }


            csptr csmat = NULL;

            csmat = malloc(sizeof(*csmat));

            coo2csptr(n, nnz, a, ia, ja, csmat);


            free(ia);
            free(a);
            free(ja);

            vbsptr vbmat = NULL;

            printf("prm->eps = %f\n",prm->eps);

            if (prm->cosine)
//                ierr = pablo( csmat, prm->eps, 2.0, &nB, &nBlock, &perm);
                ierr = init_blocks( csmat, &nBlock, &nB, &perm, prm->eps);//int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps)//parms_PCSetup(pc);
            else
                ierr = init_blocks_density( csmat, &nBlock, &nB, &perm, prm->eps);

            printf("nBlock value is %d.\n", nBlock);//%f %p %s %c
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

            outputvbmat(vbmat,"vbmatparms.coo",1);

            blocksize = (double)csmat->n / (double)nBlock;
            Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;
            if (myid == 0)
            {
                printf("\n Bsize=%-7f, Bdensity=%-7f\n",blocksize, Bdensity);

            }


            int nbb, *bia, *bja;
            BData *ba;
            bia = (int*)malloc((vbmat->n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
            bja = (int*)malloc(nnz*sizeof(int));
            ba = (BData*)malloc(nnz*sizeof(BData));
            vbsptr2colunms(vbmat, &nbb, bia, bja, ba);//int vbsptr2colunms(vbsptr mat, int *n, int *ia, int *ja)

            bja = (int*)realloc(bja, nbb*sizeof(int));


            tib4 =  parms_TimerGet(tm);
            printf("\ntime on csrvbsrC_new=%f, time on total process = %f\n",tib4-tib3, tib4);
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

            print_mem("end of proc 1");
            cleanCS( csmat );
            cleanVBMat( vbmat );
        }
        tib2 =  parms_TimerGet(tm);
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
        printf("nzding %d on %d\n", nzding[myid], myid);
        print_mem("after free");


        tib1 =  parms_TimerGet(tm);
        printf("\ntime on send=%f\n",tib1 - tib2);

        /*-------------------- Create map object */
        parms_MapCreateFromPtr(&map, n, dom, idom, MPI_COMM_WORLD, 1, NONINTERLACED);

        parms_Map_Assign_blockstructure(map, nB);//int parms_Map_Assign_blockstructure(parms_Map self, int *nB)//


        nloc = parms_MapGetLocalSize(map);


        /* Free dom and idom */

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
        ia = malloc(nzding[myid]*sizeof(*ia));
        if (ia == NULL) {
            fprintf(stderr, "cannot allocate memory for ia\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        ierr = mm_partial_read_mtx_crd_data_new(mtxfile, n, nc, nnz, ia, ja, a, idom, dom, perm, nB, nBlock, matcode);
        fclose(mtxfile);

        printf("n value is %d, nc value is %d\n", n, nc);//%f %p %s %c


        free(dom);
        free(idom);

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
            if ((mtxfile = fopen(currhs, "r")) == NULL) {
                fprintf(stderr, "Error opening rhs file\n");
                MPI_Finalize();
                exit(1);
            }

            ierr = mm_read_banner(mtxfile, &matcode);
            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read banner\n");
                fprintf(stderr, "filename = %s\n", currhs);
                MPI_Finalize();
                exit(1);
            }

            ierr = mm_read_mtx_array_size(mtxfile, &n, &nc);

            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read array size\n");
                fprintf(stderr, "filename = %s\n", currhs);
                MPI_Finalize();
                exit(1);
            }

            nrhs = 1;
            ierr = mm_read_array_data(mtxfile, n, rhstmp, matcode);
            fclose(mtxfile);

            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read rhs\n");
                fprintf(stderr, "filename = %s\n", currhs);
                MPI_Finalize();
                exit(1);
            }
        }
        else nrhs = 0;

        csptr csmat = NULL;

        csmat = malloc(sizeof(*csmat));

        coo2csptr(n, nzding[myid], a, ia, ja, csmat);//nzding is the local length array
        printf("nzding[myid] value is %d, myid = %d\n", nzding[myid], myid);//%f %p %s %c


        free(ia);
        free(a);
        free(ja);

        vbsptr vbmat = NULL;
        /*--------------------create a timer */
        tib1 =  parms_TimerGet(tm);

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

        FLOAT *rhstmpp;
        rhstmpp = (FLOAT*)malloc(n*sizeof(FLOAT));


        printf("n value is %d\n", n);//%f %p %s %c


        for( i = 0; i < n; i++ )
            rhstmpp[perm[i]] = rhstmp[i];

        int nbb, *bia, *bja;
        BData *ba;
        bia = (int*)malloc((vbmat->n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
        bja = (int*)malloc(nzding[myid]*sizeof(int));
        ba = (BData*)malloc(nzding[myid]*sizeof(BData));
        vbsptr2colunms(vbmat, &nbb, bia, bja, ba);//int vbsptr2colunms(vbsptr mat, int *n, int *ia, int *ja)

        bja = (int*)realloc(bja, nbb*sizeof(int));
        print_mem("end of partial load");

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
        parms_MatSetValues_b(A, nBlock, bim, bia, bja, ba, INSERT);//parms_MatSetValues(A, n, im, ia, ja, a, INSERT);


        free(ba);
        free(bja);
        free(bia);

        int llsize;

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


        printf("nrhs value is %d.\n", nrhs);//%f %p %s %c

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
            fprintf(fp, "The time for finding blocks=%-7f, The total blocking time=%-7f\n", tib1, tib4);
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
            fprintf(fout, "n = %d, nnz = %d\n", n, nnz);
            fprintf(fout, "The blocking information Bsize=%-7f, Bdensity=%-7f\n", blocksize, Bdensity);
            fprintf(fout, "The time for finding blocks=%-7f, The total blocking time=%-7f\n", tib1, tib4);
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


        /*--------------------Free memories */
        cleanCS( csmat );
        cleanVBMat( vbmat );//int cleanVBMat( vbsptr vbmat )
        free(x);
        free(y);
        free(rhs);
        free(perm);
        parms_MatFree_b(&A);
        parms_MapFree(&map);
        parms_PCFree_b(&pc);
        parms_SolverFree_b(&solver);
        parms_TimerFree(&tm);
        //printf("after cleanning");

        /*----Goto next matrix ---*/
//        label1000://for block structure detection
        mat++;

    }

    /*--------------Free prm--------------*/
    free(prm);
    /*--------------------Exit MPI environment */
    MPI_Finalize();
    return 0;
}
