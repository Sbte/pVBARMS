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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(__ICC)
#include <mathimf.h>
#else
#include <math.h> 
#endif 
#include "parms.h"
#include "aux.h"
//#include <sys/time.h>//new
//#include <time.h>//new
//#include "/home/p264298/cmm/PARMSproject/pARMS_3.2/src/DDPQ/globheads.h"//new
//#include "parms_mat_impl.h"//new





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

    /* Viewer object for solver */
    //  parms_Viewer  sv;



    //    FLOAT atest = 5.0 + 3.0 * I, btest;

    //    btest = conj(atest) ;

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
        printf("type value is %s\n", type);//%f %p %s %c
        //        printf("title value is %s\n", title);//%f %p %s %c
        printf("key value is %s\n", key);//%f %p %s %c
        //        printf("key value is %s\n", key);//%f %p %s %c

        //        int rsa;//to judge if the matrix is symm or skew-symm or Hermition
        int f_nnz;
        if( type[1] == 'S' || type[1] == 's' || type[1] == 'H' || type[1] == 'h')
            f_nnz = nnz+nnz-n;
        else if (type[1] == 'Z' || type[1] == 'z')
            f_nnz = nnz + nnz;
        else f_nnz = nnz;

        //        type[1] == 'S' || type[1] == 's' ||
        //                type[1] == 'Z' || type[1] == 'z' ? f_nnz = nnz + nnz - n: f_nnz = nnz;

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
        /*
  wreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],&len,a,ja,ia,rhstmp,&nrhs,
        guesol,&n,&nc,&nnz,title,key,type,&ierr);
*/
#if defined(DBL_CMPLX)
        zreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
                  guesol,&n,&nc,&nnz,title,key,type,&ierr);
#else
        readmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
                 guesol,&n,&nc,&nnz,title,key,type,&ierr);
#endif
        //output_dblvectorpa("rhstmpjustloaded",rhstmp,0, n);


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



        //        memcpy(ia, ib, (n+1)*sizeof(*ia));
        //        memcpy(ja, jb, nnz*sizeof(*ja));
        //        memcpy(a,  b, nnz*sizeof(*a));
        //        output_intvectorpa("ja",ja,0, nnz);
        //        output_intvectorpa("ia",ia,0, n+1);
        if((ierr = fullmatize(n, f_nnz, a, ja, ia, b, jb, ib, type))){
            fprintf(stderr, "cannot read matrix\n");
            MPI_Finalize();
            exit(1);
        }

//        output_dblvectorpa("BHBa", a, 0 , f_nnz);
//        output_intvectorpa("BHBja", ja, 0 , f_nnz);
//        output_intvectorpa("BHBia", ia, 0 , n+1);
        /*---------------Free CSR matrix ------------------*/

        free(ib);
        free(b);
        free(jb);

        csptr csmat = NULL;

        //mat = malloc(sizeof(*(csptr)));
        csmat = malloc(sizeof(*csmat));
  //      int colunms2csptr(int n, int *ia, int *ja, FLOAT *a, csptr mat)

//        output_dblvectorpa("mtxa", a, 0 , nnz);
//        output_intvectorpa("mtxja", ja, 0 , nnz);
//        output_intvectorpa("mtxia", ia, 0 , n+1);
        colunms2csptr(n, ia, ja, a, csmat);
        outputcsmat(csmat,"BHBmat.coo",1);

//        exit(1);

        free(ia);
        free(a);
        free(ja);

        /* ia = ib; */
        /* ja = jb; */
        /* a = b; */

        //   double t1;
        //   t1 = sys_timer();
        //csptr csmat = NULL;
        //ddouble aa;
        vbsptr vbmat = NULL;

        int nBlock, *nB = NULL, *perm = NULL;
        double tib1, tib2, tib3, tib4, blocksize;
        /*--------------------create a timer */
        parms_TimerCreate(&tm);
        parms_TimerReset(tm);


        printf("prm->eps = %f\n",prm->eps);
        //~ ierr = init_blocks_cliques( csmat, &nBlock, &nB, &perm, prm->eps);
        if (prm->cosine)
            ierr = init_blocks( csmat, &nBlock, &nB, &perm, prm->eps);//int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps)//parms_PCSetup(pc);
        else
            ierr = init_blocks_density( csmat, &nBlock, &nB, &perm, prm->eps);
        //printf("nBlock: %d\n", nBlock);getchar();
        tib1 =  parms_TimerGet(tm);
        printf("\ntime on init=%f\n",tib1);
        if(ierr != 0) {
            //fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "*** in init_blocks ierr != 0 ***\n");
            MPI_Finalize();
            exit(1);
        }
        //output_intvector("perm.coo",perm,0, n);getchar();
        if( dpermC( csmat, perm ) != 0 ) {
            fprintf( stderr, "*** dpermC error ***\n" );
            MPI_Finalize();
            exit(1);
        }
        tib2 =  parms_TimerGet(tm);
        printf("\ntime on dpermC=%f\n",tib2-tib1);
        /*-------------------- convert to block matrix. */
        vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
        //ierr = csrvbsrC( 1, nBlock, nB, csmat, vbmat );
        ierr = csrvbsrC_new( 1, nBlock, nB, csmat, vbmat );
        //outputvbmatpa(vbmat,"vbmat2.coo",1);//int outputvbmat( vbsptr vbmat, char *filename, int onebase)
        tib3 =  parms_TimerGet(tm);
        printf("\ntime on csrvbsrC_new=%f\n",tib3-tib2);
        outputvbmatpa(vbmat,"vbmat2.coo",1);


        //exit(1);
        //if (myid == 0)
        //outputvbmatpa(vbmat,"vbmat.coo",1);//int outputvbmat( vbsptr vbmat, char *filename, int onebase)

        //sprintf(buf, "%8.2fs", tib);
        //sprintf(buf, "%8.2fs", tpc);
        //fprintf(fp, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
        //fprintf(fp, "The time cost for init_blocks %7s %-s\n", "=", strtok(buf, " "));
        blocksize = (double)csmat->n / (double)nBlock;
        Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;
        if (myid == 0)
            printf("\n Bsize=%-7f, Bdensity=%-7f\n",blocksize, Bdensity);
        //MPI_Finalize();
        //goto label1000;//continue; // only for blocking information

        FLOAT *rhstmpp;
        rhstmpp = (FLOAT*)malloc(n*sizeof(FLOAT));

        //output_intvectorpa("perm_init.coo",perm,0, n);
        // output_dblvectorpa("rhstmp",rhstmp, 0, n);

        for( i = 0; i < n; i++ )
            rhstmpp[perm[i]] = rhstmp[i];
        //output_dblvectorpa("rhstmpp",rhstmpp, 0, n);

        //cleanCS(csmat);//int cleanCS(csptr amat)
        int nbb, *bia, *bja;
        BData *ba;
        bia = (int*)malloc((vbmat->n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
        bja = (int*)malloc(f_nnz*sizeof(int));
        ba = (BData*)malloc(f_nnz*sizeof(BData));
        vbsptr2colunms(vbmat, &nbb, bia, bja, ba);//int vbsptr2colunms(vbsptr mat, int *n, int *ia, int *ja)
        //iia = (int*)malloc(n*sizeof(int)   printf("\n Bsize=%-7f\n",blocksize););//w =(double*)malloc(n*sizeof(double));
        //jja = (int*)malloc(nnz*sizeof(int));
        //csptr2colunms(csmat, &nn, iia, jja);//int csptr2colunms(csptr mat, int *n, int *ia, int *ja)//
//        bja = (int*)realloc(bja, nbb*sizeof(int));
        //printf("bja[0]=%d",bja[0]);
        //output_intvector("bia.coo",bia,0, nBlock);getchar();
        //output_intvector("bja.coo",bja,0, nbb);getchar();


        tib4 =  parms_TimerGet(tm);
        printf("\ntime on csrvbsrC_new=%f, time on total process = %f\n",tib4-tib3, tib4);
        /*--------------------get the elapsed time spent on creating PC */

        //int dpermC(csptr mat, int *perm)


        //printf("\n npro=%d\n",npro);
        //int nn, *iia, *jja;
        //iia = (int*)malloc((n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
        //jja = (int*)malloc(nnz*sizeof(int));
        //csptr2colunms(csmat, &nn, iia, jja);//int csptr2colunms(csptr mat, int *n, int *ia, int *ja)//
        //  im = (int *)malloc(n*sizeof(int));
        //printf("nn = %d\n", nn );getchar();
        //output_intvector("iia.coo",iia,0, nn+1);getchar();
        //output_intvector("jja.coo",jja,0, nnz);getchar();




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
        /*
    for (i = 0; i < n; i++) {
      dom[i] = i+1;
    }
    for(i=1; i<npro; i++)
    idom[i] = idom[i-1]+(n/npro);
    idom[npro+1] = n+1;
printf("idom[%d] = %d; dom[%d] = %d \n", myid, idom[myid], idom[myid]-1,dom[idom[myid]-1]);
*/
        /*-------------------- Create map object */
        parms_MapCreateFromPtr(&map, n, dom, idom, MPI_COMM_WORLD, 1, NONINTERLACED);
        //map->nB = nB;
        parms_Map_Assign_blockstructure(map, nB);//int parms_Map_Assign_blockstructure(parms_Map self, int *nB)//
        //output_intvector("map->nB.coo",map->nB,0, nBlock);getchar();


        nloc = parms_MapGetLocalSize(map);

        //output_intvectorpa("mapperm_p",map->lvars,0, nloc);getchar();

        //getlocalsize
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
        parms_MatSetValues_b(A, nBlock, bim, bia, bja, ba, INSERT);//parms_MatSetValues(A, n, im, ia, ja, a, INSERT);

        //parms_MatSetValues_b(A, n, im, ia, ja, a, INSERT);//int parms_MatSetValues_b(parms_Mat self, int m, int *im, int *ia, int *ja, BData *values, INSERTMODE mode)
        //  parms_MatSetElementMatrix(A, n, im, ia, ja, a, INSERT);
        //  parms_MatAssembleElementMatrix(A);
        /*-------------------- free the matrix stored in CSR format(a,ja,ia) */
        /* free(a); */
        /* free(ja); */
        /* free(ia); */


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
        //outputcsmatpa(A->aux_data,"aux_data",1);//int outputcsmat ( csptr mat, char *filename, int onebase){
        //output_intvectorpa("localpermbefore",A->is->perm,0, nloc);getchar();


        /*-------------------- Setup the matrix and communication structure */
        parms_Mat_B_Setup(A);//parms_MatSetup(A);
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
            parms_VecSetValues_b(rhs, nBlock, bim, rhstmpp, INSERT, map);//changed//parms_VecSetValues_b(rhs, n, im, rhstmp, INSERT, map);
        }

        free(rhstmp);
        free(rhstmpp);
        free(im);

        free(bim);

        //output_dblvectorpa("rhs",rhs,0, llsize);getchar();

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

        //output_dblvectorpa("resvecdbl",resvec,0, llsize);getchar();
        parms_VecGetNorm2_b(resvec, &norm, map);

        //printf("norm = %20.16e\n",norm);
        free(resvec);
        /*--------------------Create preconditioner based on the matrix A. */
        parms_PCCreate_b(&pc, A);

        /*--------------------set parameters for pc */
        set_pc_params_b(pc, prm);


        /*--------------------reset the timer */
        parms_TimerReset(tm);
        parms_PCSetup_b(pc);


        /*   parms_Viewer v; */
        /*   parms_ViewerCreate(&v, "PC_Bdd");  */
        /*   parms_PCView(pc, v); */

        /*     FLOAT *xx, *yy; */
        /*     PARMS_NEWARRAY(xx, llsize); */
        /*     PARMS_NEWARRAY(yy, llsize); */
        /*     for (i=0; i < llsize; ++i) */
        /*         xx[i] = 1.0; */

        /*     /\* for(i=0;i<llsize;i++) {  *\/ */
        /*     /\*   xx[i]=(double)(rand()/(1000000));  *\/ */
        /*     /\*   //printf("xx[i] = %d \n",j);  *\/ */
        /*     /\* }  *\/ */
        /*     output_dblvectorpa("pcinitial_Bdd",xx, 0, llsize); */

        /*     parms_MatVec(A, xx, yy); */

        /*     //yy = rhs; */

        /*     for (i=0; i < llsize; ++i) */
        /*         xx[i] = 0.0; */
        /*     parms_PCApply(pc, yy, xx); */
        /*     output_dblvectorpa("pcsol_Bdd",xx, 0, llsize); */
        /*     output_dblvectorpa("pcmv_Bdd",yy, 0, llsize); */

        /* printf("after testing\n"); */
        /*   /\* MPI_Barrier(MPI_COMM_WORLD); *\/ */
        /*   /\* exit(1); *\/ */


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
        //parms_SolverSetup(solver);

        /*--------------------restart the timer */
        parms_TimerRestart(tm);

        /*--------------------Solve the linear equation */
        parms_SolverApply_b(solver, rhs, x);

        //MPI_Barrier(MPI_COMM_WORLD);
        //exit(1);

        /*--------------------get total time spent on creating the pc and solving the linear
     system */
        ttol = parms_TimerGet(tm);

        printf("The time for solving  %f in proc %d\n", ttol-tpc, myid);
        printf("The total time cost  %f in proc %d\n", ttol, myid);

        /*--------------------Get the number of iterations */
        its = parms_SolverGetIts(solver);

        /*--------------------Compute the residual error  */
        parms_MatVec(A, x, y);
        for(i=0; i<llsize; i++)//nloc
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


        //printf("after getratio\n");
        //MPI_Barrier(MPI_COMM_WORLD);
        //exit(1);
        /*
  parms_ViewerCreate(&sv, "solver.out");
  parms_TimerView(tm, sv);
*/
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
        /*
  parms_ViewerFree(&sv);
*/
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
