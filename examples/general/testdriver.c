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
#include "mmio.h"

typedef struct parms_dvcsr {
    /*! \var diag_mat The diagonal matrix stored in vcsr format.
   */
    parms_vcsr    diag_mat;
    /*!  \var offd_mat The off-diagonal matrix stored in vcsr format.
  */
    parms_vcsr    offd_mat;
    /*! \var mvhandler The parms_Comm object for the matrix-vector product.
   */
    parms_bvcsr   b_diag_mat;//vbsptr
    parms_bvcsr   b_offd_mat;//vbsptr
    parms_Comm    mvhandler;

} *parms_dvcsr;

int main(int argc, char *argv[])
{

    /* declarations related to Harwell-boeing format for reading the HB
     matri. Second part is related to I/O parameters */
    char mname[MAX_MAT][MAX_LINE],  key[MAX_LINE];//guesol[2], title[72],, type[3]
    int nrhs, nc, n, nnz, mat;// job, tmp0, tmp, tmp2, tmp3,
    int myid, ierr, i, nloc;
    /* memory usage of the preconditioning matrix */
    double ratio;
    /* working array for reading matrix */
    double norm, res1, tpc, ttol;
    FLOAT *a, *rhstmp;
    int *ja, *ia;
    int npro,its, *im;
    fprm prm;
    FILE *fp=NULL, *fout=NULL, *mtxfile=NULL;
    char *name, *iluname, *curmat, *currhs, *curname = NULL, buf[40];

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

    /*-------------------- initialize MPI environment */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &npro);

//    tmp0 = 0;
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

        if (myid == 0)
            printf("Reading matrix %s\n", curmat);

        //~ #if defined(DBL_CMPLX)
        //~ zreadmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
        //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
        //~ #else
        //~ readmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
        //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
        //~ #endif

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
        rhstmp = malloc(n*sizeof(*rhstmp));
        if (rhstmp == NULL) {
            fprintf(stderr, "cannot allocate memory for rhstmp\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
//        if(nrhs != 0)
//            tmp = 3;
//        else
//            tmp = 2;
//        tmp2 = n;
//        tmp3 = nnz;

        /*-------------------- Array sizes determined. Now call
                         wreadmtc again for really reading */
        /*
  wreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],&len,a,ja,ia,rhstmp,&nrhs,
        guesol,&n,&nc,&nnz,title,key,type,&ierr);
*/
        //~ #if defined(DBL_CMPLX)
        //~ zreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
        //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
        //~ #else
        //~ readmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
        //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
        //~ #endif

        ierr = mm_read_mtx_crd_data(mtxfile, n, nc, nnz, ia, ja, a, matcode);
        fclose(mtxfile);
//        output_intvectorpa("MTXia", ia, 0 , n+1);

        if(ierr != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "cannot read matrix\n");
            fprintf(stderr, "filename = %s\n", mname[mat]);
            MPI_Finalize();
            exit(1);
        }

        currhs = strtok(NULL, " ");
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

        /*-------------Convert from COO to CSR format ------------*/
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
//        job = 1;
#if defined(DBL_CMPLX)
        zcoocsr_(&n, &nnz, a, ia, ja, b, jb, ib);
#else
        coocsr_(&n, &nnz, a, ia, ja, b, jb, ib);
#endif
        /*---------------Free COO matrix ------------------*/

        free(ia);
        free(a);
        free(ja);

        ia = ib;
        ja = jb;
        a = b;


        output_dblvector("MTXa", a, 0 , nnz);
        output_intvector("MTXja", ja, 0 , nnz);
        output_intvector("MTXia", ia, 0 , n+1);

        output_csrmatrix("csrmat", ia, ja, a, n);
        //        output_csrmatrix("csrmat", ia, ja, a, n);

//        exit(1);
        //printf("\n npro=%d\n",npro);

//return 0;
        idom = malloc((npro+1)*sizeof(*idom));
        if (idom == NULL) {
            fprintf(stderr, "cannot allocate memory for idom\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        dom = malloc(n*sizeof(*dom));
        if (dom == NULL) {
            fprintf(stderr, "cannot allocate memory for dom\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        if (npro == 1) {
            for (i = 0; i < n; i++) {
                dom[i] = i+1;
            }
            idom[0] = 1;
            idom[1] = n+1;
        }
        else {
            riord = malloc(n*sizeof(*riord));
            if (riord == NULL) {
                fprintf(stderr, "cannot allocate memory for riord\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            mask = malloc(n*sizeof(*mask));
            if (mask == NULL) {
                fprintf(stderr, "cannot allocate memory for mask\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            jwk = malloc(2*n*sizeof(*jwk));
            if (jwk == NULL) {
                fprintf(stderr, "cannot allocate memory for jwk\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            link = malloc(n*sizeof(*link));
            if (link == NULL) {
                fprintf(stderr, "cannot allocate memory for link\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            dse_(&n, ja, ia, &npro, riord, dom, idom, mask, jwk, link);
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

        nloc = parms_MapGetLocalSize(map);

        /* Free dom and idom */
        free(dom);
        free(idom);

        /*-------------------- Create distributed vectors based on map */
        x = (FLOAT *)malloc(nloc*sizeof(FLOAT));
        rhs = (FLOAT *)malloc(nloc*sizeof(FLOAT));
        y = (FLOAT *)malloc(nloc*sizeof(FLOAT));
        resvec = (FLOAT *)malloc(nloc*sizeof(FLOAT));

        /*-------------------- create a distributed matrix based on map */
        parms_MatCreate(&A, map);
        /*-------------------- Insert values into A */
        /* Initialize global indices array */
        im = (int *)malloc(n*sizeof(int));
        for (i = 0; i < n; i++) {
            im[i] = i+1;
        }
        parms_MatSetValues(A, n, im, ia, ja, a, INSERT);
        //  parms_MatSetElementMatrix(A, n, im, ia, ja, a, INSERT);
        //  parms_MatAssembleElementMatrix(A);
        /*-------------------- free the matrix stored in CSR format(a,ja,ia) */




        free(a);
        free(ja);
        free(ia);

        /*-------------------- Setup the matrix and communication structure */
        parms_MatSetup(A);
        //parms_ViewerCreate(&v, "foomapaftersetup");
        //parms_MapView(map, v);

        //printf("A->aux_data->pj[66][0]=%d,[1]=%d \n",A->aux_data->pj[66][0],A->aux_data->pj[66][1]);
        /*-------------------- Copy rhs or Setup artifical right-hand-side vector */
        //nrhs = 0;

        //outputcsmatpa(A->aux_data,"aux_data_aftersetup",1);//int outputcsmat ( csptr mat, char *filename, int onebase){

        /* parms_Viewer v; */
        /* parms_dvcsr data; */
        /* data    = (parms_dvcsr)A->data; */
        /* parms_ViewerCreate(&v, "foomvhandler"); */
        /* parms_CommView(data->mvhandler, v);//int parms_CommView(parms_Comm self, parms_Viewer v) */
        /* //parms_ViewerFree(&v); */


        /* /\* parms_Viewer v; *\/ */
        /* parms_ViewerCreate(&v, "foomapaftersetup"); */
        /* parms_MapView(map, v); */
        /* parms_ViewerFree(&v); */

        /* MPI_Barrier(MPI_COMM_WORLD); */
        /* exit(1); */


        if(!nrhs)
        {
            for(i=0; i<nloc; i++)
            {
                x[i] = 1.0;
            }
            parms_MatVec(A, x, rhs);
        }
        else
        {
            parms_VecSetValues(rhs, n, im, rhstmp, INSERT, map);
        }

        free(rhstmp);
        free(im);

        printf("\n nloc=%d\n",nloc);
        /*--------------------Setup initial guess to a vector of all-0. */
        for(i=0; i<nloc; i++)
        {
            x[i] = 0.0;
        }
        /*--------------------Get 2-norm of initial residual */
        parms_MatVec(A, x, resvec);

        for(i=0; i<nloc; i++)
        {
            resvec[i] = resvec[i] - rhs[i];
        }
        parms_VecGetNorm2(resvec, &norm, map);
        free(resvec);


        /*--------------------Create preconditioner based on the matrix A. */
        parms_PCCreate(&pc, A);

        /*--------------------set parameters for pc */
        set_pc_params(pc, prm);



        /*--------------------create a timer */
        parms_TimerCreate(&tm);
        /*--------------------reset the timer */
        parms_TimerReset(tm);
        parms_PCSetup(pc);


        // int j;
        //for(i=0;i<10;i++) {
        //  j=(int)(rand()/(1000000));
        // printf("j = %d \n",j);
        //}
        /*   parms_Viewer v; */
        /*   parms_ViewerCreate(&v, "PC");  */
        /*   parms_PCView(pc, v); */


        /*     FLOAT *xx, *yy; */
        /*     PARMS_NEWARRAY(xx, nloc); */
        /*     PARMS_NEWARRAY(yy, nloc); */
        /*     for (i=0; i < nloc; ++i) */
        /*         xx[i] = 1.0; */

        /*     /\* for(i=0;i<nloc;i++) { */
        /*     /\*   xx[i]=(double)(rand()/(1000000));  *\/ */
        /*     /\*   //printf("xx[i] = %d \n",j);  *\/ */
        /*     /\* }  *\/ */
        /*     output_dblvectorpa("pcinitial",xx, 0, nloc); */

        /*     parms_MatVec(A, xx, yy); */
        /*     //yy = rhs; */
        /*     for (i=0; i < nloc; ++i) */
        /*         xx[i] = 0.0; */
        /*     parms_PCApply(pc, yy, xx); */
        /*     output_dblvectorpa("pcsol",xx, 0, nloc); */
        /*     output_dblvectorpa("pcmv",yy, 0, nloc); */

        /* printf("after testing\n"); */
        /*   /\* MPI_Barrier(MPI_COMM_WORLD); *\/ */
        /*   /\* exit(1); *\/ */

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
        parms_SolverSetType(solver, SOLFGMRES);

        /*--------------------set parameters for solver */
        set_solver_params(solver, prm);

        /*--------------------set up solver -- no longer needed - DOK */
        //  parms_SolverSetup(solver);

        /*--------------------restart the timer */
        parms_TimerRestart(tm);

        /*--------------------Solve the linear equation */
        parms_SolverApply(solver, rhs, x);

        /*--------------------get total time spent on creating the pc and solving the linear
     system */
        ttol = parms_TimerGet(tm);

        printf("The time for solving  %f in proc %d\n", ttol-tpc, myid);
        printf("The total time cost  %f in proc %d\n", ttol, myid);

        /*--------------------Get the number of iterations */
        its = parms_SolverGetIts(solver);

        /*--------------------Compute the residual error  */
        parms_MatVec(A, x, y);
        for(i=0; i<nloc; i++)
        {
            y[i] = rhs[i] - y[i];
        }
        parms_VecGetNorm2(y, &res1, map);
        /*--------------------processor 0 outputs the result */
        if (myid == 0) {
            parms_PCGetName(pc, &name);
            parms_PCILUGetName(pc, &iluname);

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
//            if (iluname == "VBARMS")
            if (strcmp (iluname, "VBARMS"))
                fprintf(fout, "It is in local graph compression mode\n");
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

        /*
  parms_ViewerCreate(&sv, "solver.out");
  parms_TimerView(tm, sv);
*/
        /*--------------------Free memories */
        free(x);
        free(y);
        free(rhs);
        parms_MatFree(&A);
        parms_MapFree(&map);
        parms_PCFree(&pc);
        parms_SolverFree(&solver);
        parms_TimerFree(&tm);

        /*
  parms_ViewerFree(&sv);
*/
        /*----Goto next matrix ---*/
        mat++;

    }

    /*--------------Free prm--------------*/
    free(prm);
    /*--------------------Exit MPI environment */
    MPI_Finalize();
    return 0;
}
