/*----------------------------------------------------------------------
 *                           Program dd-bin-dse
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



int main(int argc, char *argv[])
{

    /* declarations related to Harwell-boeing format for reading the HB
     matri. Second part is related to I/O parameters */
    char mname[MAX_MAT][MAX_LINE], key[MAX_LINE];//, type[3];//guesol[2], title[72],
    int nrhs, mat;//tmp0,

    long int m;    // = header(1); //row;
    long int n;    // = header(2); //column;
    long int nnz;    // = header(3);%nnz
    long int header; // header to store file type, mat or rhs;
    long int numread;
    int n_short;
    int myid, i, nloc;
    /* memory usage of the preconditioning matrix */
    double ratio;
    /* working array for reading matrix */
    double norm, res1, tpc, ttol;
    FLOAT *a, *rhstmp;
    int *ja, *nnzptr, *ia, *gia, *nzding = NULL;
    int npro,its, *im;//bim is for block case
    fprm prm;
    FILE *fp=NULL, *fout=NULL, *binfile=NULL;
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
    MPI_Init(&argc, &argv);//initilization
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);//return the pid
    MPI_Comm_size(MPI_COMM_WORLD, &npro);//return the number of the processors

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

    /* variable "mname" stores the name of the file in HB format. Read a
     Harwell-Boeing matrix. using wreadmtc c-version of sparsit
     routine - call wreadmtc a first time to determine sizes of
     arrys. read in values on the second call.
  */

    /* --- Begin loop over matrices ----*/
    mat = 0;
    while(mname[mat][1] != '#'){
        a = NULL; ja = NULL; ia = NULL; nnzptr = NULL; rhstmp = NULL;
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

        gia = malloc((m+1)*sizeof(*gia)); //nnz of each row
        gia[0] = 0;

        long int sum_nz = 0;
        for (i = 0; i < m; ++i){
            nnzptr[i] = bswap_32(nnzptr[i]);
            gia[i+1] = gia[i] + nnzptr[i];
            sum_nz += nnzptr[i];
        }
    currhs = strtok(NULL, " ");

    if (myid == 0)
    {
        ia = malloc((m+1)*sizeof(*ia)); //nnz of each row
        ia[0] = 1;
        for (i = 0; i < m; i++)
            ia[i+1] = nnzptr[i] + ia[i];

        if(sum_nz != nnz || gia[m] != nnz){
            fprintf(stderr," No-Nonzeros sum-rowlengths do not match %ld %ld %d\n", nnz, sum_nz, gia[m]);
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
            ja[i] = bswap_32(ja[i]) + 1;//swap between big endian and small endian


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
            a[i] = byteswap_double(a[i]);//need to be optimized

        /*-------------------- Array sizes determined. Now call
                         wreadmtc again for really reading */

        //  ierr = mm_read_bin_crd_data(binfile, n, nc, nnz, ia, ja, a, matcode);
        //~ fseek(binfile, 0, 0);

        fclose(binfile);
        if ((binfile = fopen(curmat, "r")) == NULL) {
            fprintf(stderr, "Error opening matrix file\n");
            fprintf(stderr, "filename = %s\n", curmat);
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
            fprintf(fp, "n = %ld, nnz = %ld\n", n, nnz);
        }


        /* output_intvector("ia.coo",ia,0, nnz); */
        /* output_intvector("ja.coo",ja,0, nnz); */
        /* output_dblvector("a.coo",a,0, nnz); */



        /* #if defined(DBL_CMPLX) */
        /*     zcoocsr_(&n, &nnz, a, ia, ja, b, jb, ib); */
        /* #else */
        /*     coocsr_(&n, &nnz, a, ia, ja, b, jb, ib); */
        /* #endif */
        /*---------------Free COO matrix ------------------*/
        /* /\*-------------Convert from COO to C-CSR format ------------*\/ */

        //        csptr csmat = NULL;

        //mat = malloc(sizeof(*(csptr)));
        //        csmat = malloc(sizeof(*csmat));

        //    output_intvectorpa("nnzptr.coo",nnzptr,0, n);
        //    exit(1);

        //        bincols2csptr(n, nnzptr, ja, a, csmat);
        //    int bincols2csptr(int n, int *nnzptr, int *ja, FLOAT *a, csptr mat)

        //    coo2csptr(n, nnz, a, ia, ja, csmat);
        //    outputcsmat(csmat,"csmat.coo",1);
        //    output_dblvectorpa("rhstmp",rhstmp,0, n);

        //    exit(1);
        //        ia = malloc(m*sizeof(*ia)); //nnz of each row

        //        output_intvector("ia.coo",ia,0, n+1);
        //        output_intvector("ja.coo",ja,0, nnz);
        //        output_dblvector("a.coo",a,0, nnz);
        //        exit(1);
        //        free(a);
        //        free(ja);

        /* ia = ib; here*/
        /* ja = jb; */
        /* a = b; */

        //        free(ia);
        //        free(a);
        //        free(ja);

        //        ia = ib;
        //        ja = jb;
        //        a = b;

        //printf("\n npro=%d\n",npro);

        n_short = (int)n;
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
            dse_(&n_short, ja, ia, &npro, riord, dom, idom, mask, jwk, link);
            free(riord);
            free(mask);
            free(jwk);
            free(link);
        }

        nzding = (int *) calloc(npro, sizeof(int));
        int i1, i3;
        int i2 = 0;
        for (i1 = 0; i1 < n; ++i1)
        {
            if (i1 == idom[i2+1]-1)
                i2++;
            i3 = dom[i1]-1;
            nzding[i2] += nnzptr[i3];
        }

        free(ja);
        free(a);
        free(ia);
    }

        if (myid != 0)
        {
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
            nzding = malloc(npro*sizeof(*nzding));
            if (dom == NULL) {
                fprintf(stderr, "cannot allocate memory for nzding\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
        }
        MPI_Bcast(idom, npro+1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(dom, n, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(nzding, npro, MPI_INT, 0, MPI_COMM_WORLD);

        ja = malloc(nzding[myid]*sizeof(*ja)); //column indeces of all non-zero entries
        if (ja== NULL) {
            fprintf(stderr, "cannot allocate memory for ja\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        a = malloc(nzding[myid]*sizeof(*a)); //values of all non-zero entries
        if (a == NULL) {
            fprintf(stderr, "cannot allocate memory for a\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        ia = malloc((n+1)*sizeof(*ia)); //values of all non-zero entries
        if (ia == NULL) {
            fprintf(stderr, "cannot allocate memory for ia\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        local_read_bin_data(binfile, m, n, nnz, ia, ja, a, idom, dom, gia);
        free(gia);

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

        printf("nloc %d n %ld\n", nloc, n);

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
        free(nnzptr);

        /*-------------------- Setup the matrix and communication structure */
        parms_MatSetup(A);
        //parms_ViewerCreate(&v, "foomapaftersetup");
        //parms_MapView(map, v);

        //printf("A->aux_data->pj[66][0]=%d,[1]=%d \n",A->aux_data->pj[66][0],A->aux_data->pj[66][1]);
        /*-------------------- Copy rhs or Setup artifical right-hand-side vector */
        //nrhs = 0;

        //outputcsmatpa(A->aux_data,"aux_data_aftersetup",1);//int outputcsmat ( csptr mat, char *filename, int onebase){

        //            outputcsmat(A->aux_data,"aux_data.coo",1);

        //                output_dblvectorpa("rhstmvpdd",rhstmp,0, n);

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

        rhstmp = malloc(n*sizeof(*rhstmp));
        if (rhstmp == NULL) {
            fprintf(stderr, "cannot allocate memory for rhstmp\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

        if (currhs != NULL)
        {
            if ((binfile = fopen(currhs, "r")) == NULL) {
                fprintf(stderr, "Error opening rhs file %s\n", currhs);
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
            fprintf(fout, "n = %ld, nnz = %ld\n", n, nnz);
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
        free(nzding);

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
