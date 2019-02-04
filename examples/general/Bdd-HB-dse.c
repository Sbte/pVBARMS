/*----------------------------------------------------------------------
 *                           Program Bdd-HB-dse
 *----------------------------------------------------------------------
 *
 *  In this test program, one processor reads the whole matrix
 *  from file and broadcast to other processors. The matrix is
 *  assumed to be in HB format. Matrix graph is then  partitioned
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

    /* declarations related to Harwell-boeing format for reading the HB
     matri. Second part is related to I/O parameters */
    char mname[MAX_MAT][MAX_LINE], guesol[2], title[72], key[8], type[3];
    int nrhs, nc, n, nnz, tmp0, tmp, tmp2, tmp3, job, mat;
    int myid, ierr, i, nloc;
    /* memory usage of the preconditioning matrix */
    double ratio, Bdensity;
    /* working array for reading matrix */
    double norm, res1, tpc, ttol;
    FLOAT *a, *rhstmp;
    int *ja, *ia;
    int npro,its, *im, *bim;//bim is for block case
    fprm prm;
    FILE *fp=NULL, *fout=NULL;
    char *name, *iluname, buf[40];

    /*-------------------- variables related to dse partitioning */
    int *riord, *dom, *idom, *mask, *jwk, *link;

    /*-------------------- variables related to pARMS */
    parms_Map       map;
    FLOAT       *x, *y, *rhs, *resvec;
    parms_Mat       A;
    parms_PC        pc;
    parms_Solver    solver;
    parms_Timer     tm;

    /*-------------------- initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &npro);

    tmp0 = 0;
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

    /* variable "mname" stores the name of the file in HB format. Read a
     Harwell-Boeing matrix. using wreadmtc c-version of sparsit
     routine - call wreadmtc a first time to determine sizes of
     arrys. read in values on the second call.
  */

    /* --- Begin loop over matrices ----*/
    mat = 0;
    while(mname[mat][1] != '#'){
        a = NULL; ja = NULL; ia = NULL; rhstmp = NULL;

#if defined(DBL_CMPLX)
        zreadmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
                  guesol,&n,&nc,&nnz,title,key,type,&ierr);
#else
        readmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
                 guesol,&n,&nc,&nnz,title,key,type,&ierr);
#endif
        int f_nnz;
        if( type[1] == 'S' || type[1] == 's' || type[1] == 'H' || type[1] == 'h')
            f_nnz = nnz+nnz-n;
        else if (type[1] == 'Z' || type[1] == 'z')
            f_nnz = nnz + nnz;
        else f_nnz = nnz;


        a = malloc(f_nnz*sizeof(*a));
        if (a == NULL) {
            fprintf(stderr, "cannot allocate memory for a\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        ja = malloc(f_nnz*sizeof(*ja));
        if (ja == NULL) {
            fprintf(stderr, "cannot allocate memory for ja\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        ia = malloc((n+1)*sizeof(*ia));
        if (ia == NULL) {
            fprintf(stderr, "cannot allocate memory for ia\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        rhstmp = malloc(n*sizeof(*rhstmp));
        if (rhstmp == NULL) {
            fprintf(stderr, "cannot allocate memory for rhstmp\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        if(nrhs != 0)
            tmp = 3;
        else
            tmp = 2;
        tmp2 = n;
        tmp3 = nnz;

        /*-------------------- Array sizes determined. Now call
                         wreadmtc again for really reading */
#if defined(DBL_CMPLX)
        zreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
                  guesol,&n,&nc,&nnz,title,key,type,&ierr);
#else
        readmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
                 guesol,&n,&nc,&nnz,title,key,type,&ierr);
#endif


        if(ierr != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "cannot read matrix\n");
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
            fprintf(fp, "\nMatrix %d: %.*s %.*s \n",(mat+1),8,key,3,type);
            fprintf(fp, "n = %d, nnz = %d\n", n, nnz);
        }



        /*-------------Convert from CSC to CSR format ------------*/
        int *jb, *ib;
        FLOAT *b;
        b   = malloc(nnz*sizeof(*b));
        if (b == NULL) {
            fprintf(stderr, "cannot allocate memory for b\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        jb  = malloc(nnz*sizeof(*jb));
        if (jb == NULL) {
            fprintf(stderr, "cannot allocate memory for jb\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        ib  = malloc((n+1)*sizeof(*ib));
        if (ib == NULL) {
            fprintf(stderr, "cannot allocate memory for ib\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        job = 1;
#if defined(DBL_CMPLX)
        zcsrcsc_(&n, &job, &job, a, ja, ia, b, jb, ib);
#else
        csrcsc_(&n, &job, &job, a, ja, ia, b, jb, ib);
#endif
        /*---------------Copy CSR matrix ------------------*/



        if((ierr = fullmatize(n, f_nnz, a, ja, ia, b, jb, ib, type))){
            fprintf(stderr, "cannot read matrix\n");
            MPI_Finalize();
            exit(1);
        }

        /*---------------Free CSR matrix ------------------*/

        free(ib);
        free(b);
        free(jb);

        csptr csmat = NULL;

        csmat = malloc(sizeof(*csmat));

        colunms2csptr(n, ia, ja, a, csmat);
        outputcsmat(csmat,"BHBmat.coo",1);


        free(ia);
        free(a);
        free(ja);

        vbsptr vbmat = NULL;

        int nBlock, *nB = NULL, *perm = NULL;
        double tib1, tib2, tib3, tib4, blocksize;
        /*--------------------create a timer */
        parms_TimerCreate(&tm);
        parms_TimerReset(tm);


        printf("prm->eps = %f\n",prm->eps);
        if (prm->constant_block_size > 0)
            ierr = init_blocks_constant( csmat, &nBlock, &nB, &perm, prm->constant_block_size);
        else if (prm->cosine)
            ierr = init_blocks( csmat, &nBlock, &nB, &perm, prm->eps);
        else
            ierr = init_blocks_density( csmat, &nBlock, &nB, &perm, prm->eps);
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
        outputvbmatpa(vbmat,"vbmat2.coo",1);


        blocksize = (double)csmat->n / (double)nBlock;
        Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;
        if (myid == 0)
            printf("\n Bsize=%-7f, Bdensity=%-7f\n",blocksize, Bdensity);

        FLOAT *rhstmpp;
        rhstmpp = (FLOAT*)malloc(n*sizeof(FLOAT));

        for( i = 0; i < n; i++ )
            rhstmpp[perm[i]] = rhstmp[i];

        int nbb, *bia, *bja;
        BData *ba;
        bia = (int*)malloc((vbmat->n+1)*sizeof(int));
        bja = (int*)malloc(f_nnz*sizeof(int));
        ba = (BData*)malloc(f_nnz*sizeof(BData));
        vbsptr2colunms(vbmat, &nbb, bia, bja, ba);


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
        /*-------------------- Create map object */
        parms_MapCreateFromPtr(&map, n, dom, idom, MPI_COMM_WORLD, 1, NONINTERLACED);
        parms_Map_Assign_blockstructure(map, nB);

        nloc = parms_MapGetLocalSize(map);

        /* Free dom and idom */
        free(dom);
        free(idom);

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

        llsize = parms_MapGetLocallSize(map);
//        printf("llsize=%d, nloc = %d\n", llsize, nloc);
        /*-------------------- Create distributed vectors based on map */
        x = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        rhs = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        y = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        resvec = (FLOAT *)malloc(llsize*sizeof(FLOAT));


        /*-------------------- Setup the matrix and communication structure */
        parms_Mat_B_Setup(A);
        printf("\n npro=%d\n",npro);

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
        parms_PCGetRatio(pc, &ratio);
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
            fprintf(fout, "\nMatrix %d: %.*s %.*s \n",(mat+1),8,key,3,type);
            fprintf(fout, "n = %d, nnz = %d\n", n, nnz);
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
        /*----Goto next matrix ---*/
        //label1000://for block structure detection
        mat++;

    }

    /*--------------Free prm--------------*/
    free(prm);
    /*--------------------Exit MPI environment */
    MPI_Finalize();
    return 0;
}
