/*----------------------------------------------------------------------
 * Parallel Multi-Level Block ILUT PRECONDITIONER
 * vbarms2    :  parallel vbarms2
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#define  PERMTOL  0.99   /*  0 --> no permutation 0.01 to 0.1 good  */
#include "protos.h" 

static int vbarms_free_vcsr(parms_Operator *self)
{
    parms_barms_data data;

    data = (parms_barms_data)(*self)->data;
    cleanBARMS(data);
    //printf("data = %p\n", data);
    PARMS_FREE(data);
    return 0;
}

static int vbarms_view_vcsr(parms_Operator self, parms_Viewer v)
{
    /* So far, only available for viewing last level operator */

    parms_barms_data data;
    vbiluptr LU;
    int i, j, n, nnz, *ja, dimR, dimC, *bsz, col, k;
    BData *ba;
    FILE *fp;

    parms_ViewerGetFP(v, &fp);
    data = (parms_barms_data)self->data;
    LU = data->ilus->lu;
    bsz = data->ilus->lu->bsz;
    n = LU->n;

    /* L part */
    fprintf(fp, "L part of the last level matrix:\n");
    fprintf(fp, "n = %d\n", n);
    for (i = 0; i < n; i++) {
        dimR = B_DIM(bsz,i);
        nnz = LU->L->nzcount[i];
        ja = LU->L->ja[i];
        ba = LU->L->ba[i];
        fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
        for (j = 0; j < nnz; j++)
        {
            col = ja[j];
            dimC = B_DIM(bsz, col-LU->L->lsize);
            for (k = 0; k < dimC*dimR; k++)
            {
#if defined(DBL_CMPLX)
                fprintf(fp, "%d,%d,%f  %f\n", bsz[i]+k % dimR+1, bsz[col-LU->L->lsize]+k/dimR+1, creal(ba[j][k]), cimag(ba[j][k]));
#else
                fprintf(fp, "%d,%d,%f\n", bsz[i]+k % dimR+1, bsz[col-LU->L->lsize]+k/dimR+1, ba[j][k]);
#endif
            }
        }
    }

    /* U part */
    fprintf(fp, "U part of the matrix:\n");
    fprintf(fp, "n = %d\n", n);
    for (i = 0; i < n; i++) {
        dimR = B_DIM(bsz,i);
        nnz = LU->U->nzcount[i];
        ja = LU->U->ja[i];
        ba = LU->U->ba[i];
        fprintf(fp, "nnzrow[%d]=%d\n", i, nnz);
        for (j = 0; j < nnz; j++)
        {
            col = ja[j];
            dimC = B_DIM(bsz, col-LU->U->lsize);
            for (k = 0; k < dimC*dimR; k++)
            {
#if defined(DBL_CMPLX)
                fprintf(fp, "%d,%d,%f  %f\n", bsz[i]+k % dimR+1, bsz[col-LU->U->lsize]+k/dimR+1, creal(ba[j][k]), cimag(ba[j][k]));
#else
                fprintf(fp, "%d,%d,%f\n", bsz[i]+k % dimR+1, bsz[col-LU->U->lsize]+k/dimR+1, ba[j][k]);
#endif
            }
        }
    }
    parms_ViewerStoreFP(v, fp);

    return 0;
}

static void parms_barms_nnz(parms_Operator self, int *nnz_mat, int
                            *nnz_pc)
{
    parms_barms_data data;
//    printf("in parms_barms_nnz\n");
    data = (parms_barms_data)self->data;
    *nnz_mat = data->nnz_mat;
    *nnz_pc  = data->nnz_prec;
}

static int parms_barms_lsol_vcsr(parms_Operator self, FLOAT *y, FLOAT *x)
{
    parms_barms_data data;
    vbp4ptr levmat;
    vbilutptr ilus;
    int *ipar, schur_start, nlev, n;

    data = (parms_barms_data)self->data;
    levmat = data->levmat;
    ilus   = data->ilus;
    ipar   = data->ipar;
    schur_start = data->schur_start;//data->bsz[schur_start];

    if (ipar[0] == 0) {
        vbLsolp(schur_start, ilus->lu, y, x);
        return 0;
    }
    nlev = data->nlev;
    n = data->n;
    PARMS_MEMCPY(x, y, n);
    BLvsol2(x, nlev, levmat, ilus, 0) ;
    return 0;
}

static int parms_barms_sol_vcsr(parms_Operator self, FLOAT *y,
                                FLOAT *x)
{
    parms_barms_data data;
    int n, i;
    int *perm;
    FLOAT *xtemp;
    //printf("in barms_sol\n");

    data = (parms_barms_data)self->data;

    n = data->n;

    if (self->blockperm){
        perm = self->blockperm;

        // printf("in rhs and solution perm\n");

        //output_intvectorpa("perm_sol.coo",perm,0, n);
        //output_dblvectorpa("x_sol.coo",x,0, n);

        PARMS_NEWARRAY(xtemp, n);
        // output_dblvectorpa("xtempbeforepermutation.coo",y,0, n);

        //PARMS_MEMCPY(xtemp,y , n);
        for( i = 0; i < n; i++ )
            xtemp[perm[i]] = y[i];

        // output_dblvectorpa("xtempbefore.coo",xtemp,0, n);

        //PARMS_MEMCPY(xtemp,y , n);
        barmsol2(xtemp, data);
        //output_dblvectorpa("xtempafter.coo",xtemp,0, n);

        // PARMS_MEMCPY(xtemp,x , n);
        /* MPI_Barrier(MPI_COMM_WORLD); */
        /* exit(1); */

        //PARMS_MEMCPY(x,xtemp , n);

        for( i = 0; i < n; i++ )
            x[i] = xtemp[perm[i]];

        PARMS_FREE(xtemp);

    }
    else{
        PARMS_MEMCPY(x, y, n);
        //output_dblvectorpa("y.coo",y,0, n);

        barmsol2(x, data);
    }

    return 0;
}

static int parms_barms_invs_vcsr(parms_Operator self,  FLOAT *y, FLOAT
                                 *x)
{
    parms_barms_data data;
    vbilutptr ilus;
    int *ipar, schur_start, n;

    data = (parms_barms_data)self->data;
    schur_start = data->schur_start;
    ilus   = data->ilus;
    ipar   = data->ipar;

    n = ilus->bsz[ilus->n];
    if (y == NULL || x == NULL) {
        return 0;
    }
    if (ipar[0] == 0)  {
        vbinvsp(schur_start, ilus->lu, y, x);
    }
    else {
        if (x != y)
            PARMS_MEMCPY(x, y, n);
        if (ilus->lu)
            VBSchLUsol(ilus, x);
        else
            pSchLUsol(ilus, x);
    }

    return 0;
}

static int parms_barms_ascend_vcsr(parms_Operator self, FLOAT *y, FLOAT
                                   *x)
{
    parms_barms_data data;
    vbp4ptr levmat=NULL, last=NULL;
    vbilutptr ilus;
    int *ipar, schur_start, n, nloc, lenB, first;

    data = (parms_barms_data)self->data;
    levmat = data->levmat;
    ilus   = data->ilus;
    schur_start = data->schur_start;
    ipar   = data->ipar;
    n      = data->n;

    if (ipar[0] == 0) {
        vbUsolp(schur_start, ilus->lu, y, x);
        return 0;
    }

    while (levmat) {
        last = levmat;
        levmat = levmat->next;
    }
    levmat = last;

    nloc=levmat->bsz[levmat->n];
    lenB=levmat->bsz[levmat->nB];
    first = n - nloc;
    /*-------------------- last level                                 */
    first += lenB;
    /*-------------------- other levels                               */
    while (levmat) {
        nloc = levmat->bsz[levmat->n];
        first -= levmat->bsz[levmat->nB];
        if (levmat->n)
            vbascend(levmat, &x[first],&x[first]);
        /*-------------------- right scaling */
        if (levmat->D2 !=  NULL)
            vbdscale(nloc, levmat->D2, &x[first], &x[first]) ;
        levmat = levmat->prev;
    }
    return 0;
}

static int parms_barms_getssize_vcsr(parms_Operator self)
{
    parms_barms_data data;

    data = (parms_barms_data)self->data;
    return data->schur_start;
}

static struct parms_Operator_ops parms_barms_sol_vptr = {
    parms_barms_sol_vcsr,
    parms_barms_lsol_vcsr,
    parms_barms_invs_vcsr,
    parms_barms_ascend_vcsr,
    parms_barms_getssize_vcsr,
    parms_barms_nnz,
    vbarms_free_vcsr,
    vbarms_view_vcsr
};


int parms_barms_vcsr(parms_Mat self, parms_FactParam param, void *mat,
                     parms_Operator *op)
{
    /*---------------------------------------------------------------------
| MULTI-LEVEL BLOCK ILUT PRECONDITIONER.
| ealier  version:  June 23, 1999  BJS -- 
| version2: Dec. 07th, 2000, YS  [reorganized ]
| version 3 (latest) Aug. 2005.  [reorganized + includes ddpq]
+---------------------------------------------------------------------- 
| ON ENTRY:
| ========= 
| ( Amat ) = original matrix A stored in C-style Compressed Sparse
|            Row (CSR) format -- 
|            see LIB/heads.h for the formal definition of this format.
|
| ipar[0:17]  = integer array to store parameters for 
|       vbarms construction (arms2) 
|
|       ipar[0]:=nlev.  number of levels (reduction processes). 
|                       see also "on return" below. 
| 
|       ipar[1]:= level-reordering option to be used.  
|                 if ipar[1]==0 ARMS uses a block independent set ordering
|                  with a sort of reverse cutill Mc Kee ordering to build 
|                  the blocks. This yields a symmetric ordering. 
|                  in this case the reordering routine called is indsetC
|                 if ipar[1] == 1, then a nonsymmetric ordering is used.
|                  In this case, the B matrix is constructed to be as
|                  diagonally dominant as possible and as sparse as possble.
|                  in this case the reordering routine called is ddPQ.
|                 
|       ipar[2]:=bsize. for indset  Dimension of the blocks. 
|                  bsize is only a target block size. The size of 
|                  each block can vary and is >= bsize. 
|                  for ddPQ: this is only the smallest size of the 
|                  last level. vbarms will stop when either the number 
|                  of levels reaches nlev (ipar[0]) or the size of the
|                  next level (C block) is less than bsize.
|
|       ipar[3]:=iout   if (iout > 0) statistics on the run are 
|                       printed to FILE *ft
|
|	ipar[4]:= Krylov subspace dimension for last level 
|		    ipar[4] == 0 means only backward/forward solve
|		    is performed on last level.                  
|	ipar[5]:=  maximum # iterations on last level
|
|       ipar[6-9] NOT used [reserved for later use] - set to zero.
| 
| The following set method options for vbarms2. Their default values can
| all be set to zero if desired. 
|
|       ipar[10-13] == meth[0:3] = method flags for interlevel blocks
|       ipar[14-17] == meth[0:3] = method flags for last level block - 
|       with the following meaning in both cases:
|            meth[0] nonsummetric permutations of  1: yes. affects rperm
|                    USED FOR LAST SCHUR COMPLEMENT 
|            meth[1] permutations of columns 0:no 1: yes. So far this is
|                    USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. 
|                    (so ipar[11] does no matter - enter zero). If 
|                    ipar[15] is one then ILUTP will be used instead 
|                    of ILUT. Permutation data stored in: perm2. 
|            meth[2] diag. row scaling. 0:no 1:yes. Data: D1
|            meth[3] diag. column scaling. 0:no 1:yes. Data: D2
|       all transformations related to parametres in meth[*] (permutation, 
|       scaling,..) are applied to the matrix before processing it 
| 
| ft       =  file for printing statistics on run
|
| droptol  = Threshold parameters for dropping elements in ILU 
|            factorization.
|            droptol[0:4] = related to the multilevel  block factorization
|            droptol[5:6] = related to ILU factorization of last block.
|            This flexibility is more than is really needed. one can use
|            a single parameter for all. it is preferable to use one value
|            for droptol[0:4] and another (smaller) for droptol[5:6]
|            droptol[0] = threshold for dropping  in L [B]. See piluNEW.c:
|            droptol[1] = threshold for dropping  in U [B].
|            droptol[2] = threshold for dropping  in L^{-1} F 
|            droptol[3] = threshold for dropping  in E U^{-1} 
|            droptol[4] = threshold for dropping  in Schur complement
|            droptol[5] = threshold for dropping  in L in last block
|              [see ilutpC.c]
|            droptol[6] = threshold for dropping  in U in last block
|              [see ilutpC.c]
|             This provides a rich selection - though in practice only 4
|             parameters are needed [which can be set to be the same 
              actually] -- indeed it makes sense to take
|             droptol[0] = droptol[1],  droptol[2] = droptol[3], 
|             and droptol[4] = droptol[5]
|
| lfil     = lfil[0:6] is an array containing the fill-in parameters.
|            similar explanations as above, namely:
|            lfil[0] = amount of fill-in kept  in L [B]. 
|            lfil[1] = amount of fill-in kept  in U [B].
|            lfil[2] = amount of fill-in kept  in E L\inv 
|            lfil[3] = amount of fill-in kept  in U \inv F
|            lfil[4] = amount of fill-in kept  in S    .
|            lfil[5] = amount of fill-in kept  in L_S  .
|            lfil[6] = amount of fill-in kept  in U_S 
|             
| tolind   = tolerance parameter used by the indset function. 
|            a row is not accepted into the independent set if 
|            the *relative* diagonal tolerance is below tolind.
|            see indset function for details. Good values are 
|            between 0.05 and 0.5 -- larger values tend to be better
|            for harder problems.
| 
| ON RETURN:
|=============
|
| (PreMat)  = vbarms data structure which consists of two parts:
|             levmat and ilsch. 
|
| ++(levmat)= permuted and sorted matrices for each level of the block 
|             factorization stored in PerMat4 struct. Specifically
|             each level in levmat contains the 4 matrices in:
|
|
|            |\         |       |
|            |  \   U   |       |
|            |    \     |   F   |
|            |  L   \   |       |
|            |        \ |       |
|            |----------+-------|
|            |          |       |
|            |    E     |       |
|            |          |       |
|            
|            plus a few other things. See LIB/heads.h for details.
|
| ++(ilsch) = IluSpar struct. If the block of the last level is:
|
|                        |  B    F |
|                  A_l = |         | 
|                        |  E    C |
|
|             then IluSpar contains the block C and an ILU
|             factorization (matrices L and U) for the last 
|             Schur complement [S ~ C - E inv(B) F ]
|             (modulo dropping) see LIB/heads.h for details.
|
| ipar[0]   = number of levels found (may differ from input value) 
|
+---------------------------------------------------------------------*/
    /* local matrix object */
    parms_Operator newOpt;
    parms_barms_data data;
    parms_bvcsr      Amat;
    parms_Map       is;
    /*-------------------- function  prototyping  done in LIB/protos.h    */
    /*-------------------- move above to protos.h */
    vbp4ptr levp=NULL, levc=NULL, levn=NULL, levmat=NULL;
    vbsptr schur=NULL, B=NULL, F=NULL, E=NULL, C=NULL;
    vbilutptr ilsch=NULL;
    csptr csmat2=NULL;
    /*-------------------- local variables  (initialized)   */
    double *dd1 = NULL, *dd2 = NULL;
    double *droptol, tolind;
    int nlev, bsize, iout, ierr = 0;
    int Storedd12 = 0; //for judging dd1 and dd2 should be stored or not
    int *ipar, *lfil;
    int methL[4], methS[4];
    /*--------------------  local variables  (not initialized)   */
    int nA, nB, nC, j, n, ilev, symperm, schur_start, nbnd, nn, schur_start_p;//i,
    FILE *ft;
    /*--------------------    work arrays:    */
    int *iwork, *uwork, *iworkn, *uworkn, *nBB, nset;
    /*   timer arrays:  */
    /*   double *symtime, *unstime, *factime, *tottime;*/
    /*---------------------------BEGIN ARMS-------------------------------*/
    /*   schur matrix starts being original A */

    /*-------------------- begin                                         */
    n = param->n;
    nbnd = schur_start = param->schur_start;
    is = self->is;
    if (schur_start == -1) {
        nbnd = is->schur_start;
        printf("nbnd = %d \n",nbnd);
    }



    // This statement is invalid not only here but also in arms2 ????
    if (!param->isalloc) {
        // printf("if (!param->isalloc)\n");

        parms_OperatorCreate(&newOpt);
        PARMS_MEMCPY(newOpt->ops, &parms_barms_sol_vptr, 1);
        PARMS_NEW(data);
        PARMS_NEW(data->levmat);
        levmat = data->levmat;
        //printf("data->levmat=%p\n", data->levmat);
        PARMS_NEW(data->ilus);
        ilsch =  data->ilus;
        newOpt->data = data;
        *op = newOpt;
    }
    else {
        //printf("in else\n");
        data = (*op)->data;
        //printf("nothing wrong\n");
    }
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */
    Amat = (parms_bvcsr)mat;

    /* compute the number of nonzero entries in block form mat */
    data->nnz_mat = nnzVBMat1(Amat);//memVBMat(Amat);//nnzVBMat1(Amat);

    /* int myid; */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &myid); */
    /* printf("The number of non-zero entries is %d in local block form matrix, pid = %d \n",data->nnz_mat, myid); */
    /* printf("The number of stored entries is %d in local block form matrix, pid = %d \n",memVBMat(Amat), myid); */
    /* outputvbmatpa(Amat,"localvbmat_data",1); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */

    /* retrieve data from input objects */
    ipar = param->ipar;
    lfil = param->lfil;
//    printf("param->lfil[7] value is %d\n", param->lfil[7]);//%f %p %s %c
//    printf("lfil[7] value is %d\n", lfil[7]);//%f %p %s %c

    droptol = param->droptol;
    tolind  = param->tolind;

    nlev = ipar[0];
    //printf("nlev=%d\n", nlev);
//    printf("in barms_vcsr\n");
    fflush(stdout);
    bsize = ipar[2];
    iout = ipar[3];

    if (iout > 0 ) {
        ft = stdout;
    }
    else{
        ft = NULL;
    }

    ierr = 0;
    /*---------------------------------------------------------------------
| The matrix (a,ja,ia) plays role of Schur compl. from the 0th level.
+--------------------------------------------------------------------*/

    nC = nA = n = Amat->n;
    nn = Amat->bsz[n];

    schur_start_p = Amat->bsz[is->schur_start]; //new, for the point index of schur_start, needed in schur-solver

    //schur_start_p = Amat->bsz[is->schur_start]; //new, for the point index of schur_start, needed in schur-solver

    if (bsize >= n) bsize = n-1;
    levmat->n = n; levmat->nB = 0;

    //FIXME: Not sure why this is needed here and not in arms2
    //schur = (vbsptr) Malloc(sizeof(VBSparMat), "arms2:7" );
    //vbscpy(Amat, schur);
    schur = Amat;
    levc = levmat;
    /*--------------------------------------- */
    levc->prev = levc->next = levp = NULL;
    levc->n = 0;

    iwork = (int *) Malloc(n*sizeof(int), "barms2:2.1" );
    nBB = (int *) Malloc(n*sizeof(int), "barms2:2.2" );

    memcpy(methL, &ipar[10], 4*sizeof(int));
    memcpy(methS, &ipar[14], 4*sizeof(int));
    /*---------------------------------------------------------------------
    | The preconditioner construction is divided into two parts:
    |   1st part: construct and store multi-level L and U factors;
    |   2nd part: construct the ILUT factorization for the coarsest level
    +--------------------------------------------------------------------*/
    if ( (iout > 0)  && (nlev > 0) ) {
        fprintf(ft,"  \n");
        fprintf(ft,"Level   Total Unknowns    B-block   Coarse set\n");
        fprintf(ft,"=====   ==============    =======   ==========\n");
    }
    /*---------------------------------------------------------------------
    | main loop to construct multi-level LU preconditioner. Loop is on the
    | level ilev. nA is the dimension of matrix A_l at each level.
    +--------------------------------------------------------------------*/
    for (ilev = 0; ilev < nlev; ilev++) {
        /*-------------------- new nA is old nC -- from previous level */
        nA = nC;

        if ( nA <= bsize )  goto label1000;


        /*-------------------- allocate work space                        */
        iworkn = (int *) Malloc(schur->bsz[nA]*sizeof(int), "barms2:2.5" );
        symperm = 0;    /* 0nly needed in cleanP4 */
        if (ipar[1] == 1) {
            uworkn = (int *) Malloc(schur->bsz[nA]*sizeof(int), "barms2:2.6" );
            uwork = (int *) Malloc(n*sizeof(int), "barms2:2.6" );
        }
        else{
            symperm = 1;
            uwork = iwork;
            uworkn = iworkn;
        }
        /*-------------------- SCALING*/
        dd1 = NULL;
        dd2 = NULL;
        if (methL[2]) {
            dd1 = (double *) Malloc(schur->bsz[nA]*sizeof(double), "barms2:3" );
            j=vbfroscalC(schur, dd1,1);
            if (j) printf("ERROR in vbfroscalC -  row %d  is a zero row\n",j);
        }

        if (methL[3]) {
            dd2 = (double *) Malloc(schur->bsz[nA]*sizeof(double), "barms2:4" );
            j=vbfcoscalC(schur, dd2,1);
            if (j) printf("ERROR in vbfcoscalC - column %d is a zero column\n",j);
        }
        /*--------------------independent-sets-permutation-------------------
      |  do reordering -- The matrix and its transpose are used.
      +--------------------------------------------------------------------*/
        /* if (SHIFTTOL > 0.0) shiftsD(schur,SHIFTTOL);    */
        csmat2 = (csptr) malloc(sizeof(SparMat));
        ierr = vbsrc2csr(schur, csmat2);
        if (ipar[1] == 1)
            PQperm(csmat2, uwork, iwork, &nB, tolind, nbnd) ;
        else
            //  indsetC(csmat2, bsize, iwork, &nB, tolind, nbnd) ;
            indsetC2(csmat2, bsize, iwork, &nB, tolind, nbnd, nBB, &nset) ;
        cleanCS(csmat2);
        /*---------------------------------------------------------------------
      | nB is the total number of nodes in the independent set.
      | nC : nA - nB = the size of the reduced system.
      +--------------------------------------------------------------------*/
        nC = nA - nB;
        nbnd -= nB;
        /*   if the size of B or C is zero , exit the main loop  */
        /*   printf ("  nB %d nC %d \n",nB, nC); */
        //printf("nB=%d\n", nB);
        //printf("nC=%d\n", nC);
        if ( nB == 0 || nC == 0 )  {
            Storedd12 = 1;
            free(iworkn);iworkn = NULL;
            goto label1000;
        }
        /*int myid;
MPI_Comm_rank(MPI_COMM_WORLD, &myid);
printf("myid=%d\n", ilev);
  printf("ilev=%d\n", ilev);
  printf("nA=%d\n", nA);
  printf("bsize=%d\n", bsize);*/
        /*---------------------------------------------------------------------
      | The matrix for the current level is in (schur).
      | The permutations arrays are in iwork and uwork (row).
      | The routines rpermC, cpermC permute the matrix in place.
      *-----------------------------------------------------------------------*/
        /*   DEBUG : SHOULD THIS GO BEFORE GOTO LABEL1000 ?? */
        vbrpermC(schur,uwork);
        vbcpermC(schur,iwork);
        ierr = permuten(iwork, iworkn, schur->bsz, nA);
        /*   prtC(schur, ilev) ;   print matrix - debugging */
        /*-----------------------------------------------------------------------
      | If this is the first level, the permuted matrix is stored in
      | (levc) = (levmat).  Otherwise, the next level is created in (levc).
      +--------------------------------------------------------------------*/
        if (ilev > 0) {
            /*-   delete C matrix of any level except last one (no longer needed) */
            cleanVBMat(C);
            levn = (vbp4ptr) Malloc(sizeof(VBPer4Mat), "barms2:6" );
            /* levc->prev = levp; */
            levc->next = levn;
            levp = levc;
            levc = levn;
            levc->prev = levp;
        }
        /*-------------------- p4ptr struct for current schur complement */
        B = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:7" );
        E = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:8" );
        F = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:9" );
        C = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:10" );
        vbSplit4copy(schur, nB, nC, B, F, E, C);
        ierr = setupVBP4(levc, nB, nC, F, E, schur->bsz);
        /*--------------------     copy a few pointers       ---- */
        levc->perm  = iworkn;
        levc->rperm = uworkn;
        levc->symperm = symperm;
        levc->D1=dd1;
        levc->D2=dd2;
        /*---------------------------------------------------------------------
      | a copy of the matrix (schur) has been permuted. Now perform the
      | block factorization:
      |
      | | B   F |       | L       0 |     | U  L^-1 F |
      | |       |   =   |           |  X  |           | = L x U
      | | E   C |       | E U^-1  I |     | 0    A1   |
      |
      | The factors E U^-1 and L^-1 F are discarded after the factorization.
      |
      +--------------------------------------------------------------------*/
        if (iout > 0)
            fprintf(ft,"%3d %13d %13d %10d\n", ilev+1,schur->bsz[nA],schur->bsz[nB],schur->bsz[nA]-schur->bsz[nB]);
        /*---------------------------------------------------------------------
      | PILUT constructs one level of the block ILU fact.  The permuted matrix
      | is in (levc).  The L and U factors will be stored in the p4mat struct.
      | destroy current Schur  complement - no longer needed  - and set-up new
      | one for next level...
      +--------------------------------------------------------------------*/
        cleanVBMat(schur);
        schur = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:11" );
        /*----------------------------------------------------------------------
      | calling PILU to construct this level block factorization
      | ! core dump in extreme case of empty matrices.
      +----------------------------------------------------------------------*/
        //FIXME!!!!!
        //~ if (ilutp == 1)
        //~ ierr = computschurpartitionp(levc, B, C, droptol, lfil, ipar, schur, nBB, nset);
        //~ else
        //~ printf("levc=%p\n", levc);
        //~ ierr = computschurpartition(levc, B, C, droptol, lfil, schur, nBB, nset);
        setupVBMat(schur, nC, NULL);
        if (lfil[7] == -1)
            ierr = vbiluNEW(levc, B, C, droptol, lfil, schur);
        else
            ierr = vbilukNEW(levc, B, C, lfil[7], schur);
        //~ outputvbmat(schur, "newschur", 1);
        //~ outputvbmat(C, "cvbarms", 1);
        //ierr = computschurpartition(levc, B, C, droptol, lfil, schur, nBB, nset);
        //~ outputvbmat(schur, "oldschur", 1);
        /* prtC(levc->L, ilev) ; */
        if (ierr) {
            if (ft)
                fprintf(ft," ERROR IN  PILU  -- IERR = %d\n", ierr);
            return(1);
        }
        cleanVBMat(B);
    }
    /*---------------------------------------------------------------------
    |   done with the reduction. Record the number of levels in ipar[0]
    |**********************************************************************
    +--------------------------------------------------------------------*/
label1000:
    /* printf (" nnz_Schur %d \n",cs_nnz (schur)); */
    if (iwork) free(iwork);
    if (nBB) free(nBB);
    levc->next = NULL;
    ipar[0] = ilev;
    data->nlev = ilev;
    data->n = nn;
    nC = schur->n;
    setupVBILUT(ilsch, nC, schur->bsz);

    if ( Storedd12 )  {
        ilsch->dd1 = dd1;
        ilsch->dd2 = dd2;
    }
    /*--------------------------------------------------------------------*/
    /* define C-matrix (member of ilsch) to be last C matrix */
    if (ilev > 0) ilsch->C = C;

    /*-------------------- for ilut fact of schur matrix */
    /*  SCALING  */

    ilsch->D1 = NULL;
    if (methS[2] || methS[1] == 0) {
        //printf("in scale\n");
        ilsch->D1 = (double *) Malloc(schur->bsz[nC]*sizeof(double), "barms2:iluschD1" );
        j=vbfroscalC(schur, ilsch->D1, 1);
        if (j) printf("ERROR in vbfroscalC - row %d is a zero row\n",j);
    }

    ilsch->D2  = NULL;
    if (methS[3]) {
        ilsch->D2 = (double *) Malloc(schur->bsz[nC]*sizeof(double), "barms2:iluschD1" );
        j =vbfcoscalC(schur, ilsch->D2, 1);
        if (j) printf("ERROR in vbfcoscalC - column %d is a zero column\n",j);
    }

    /*---------------------------------------------------------------------
    |     get ILUT factorization for the last reduced system.
    +--------------------------------------------------------------------*/
    //~ uwork = NULL;
    iwork = NULL;

    ilsch->rperm = iwork; //uwork;
    ilsch->perm  = iwork;

    /*   printf("  lf : %d  %d  %d  %d  %d  %d  %d  \n",lfil[0],
       lfil[1], lfil[2], lfil[3], lfil[4], lfil[5], lfil[6]) ; */

    ierr = vbilutD(schur, droptol, lfil, ilsch);
    /*---------- OPTIMIZATRION: NEED TO COMPOUND THE TWO
    RIGHT PERMUTATIONS -- CHANGES HERE AND IN
    USCHUR SOLVE ==  compound permutations */
    if (ierr) {
        if (ft)
            fprintf(ft," ERROR IN  ILUT -- IERR = %d\n", ierr);
        return(1);
    }
    /* Last Schur complement no longer needed */
    cleanVBMat(schur);
    data->nnz_prec = nnz_vbarms(data, stderr);
    printf("bsz %d %d\n", ilsch->bsz[ilsch->n], ilsch->lu->bsz[ilsch->lu->n]);
    //~ data->ind = n - ilsch->n;
    data->ind = nn - ilsch->bsz[ilsch->n];
    //  printf("n=%d, ilschn=%d\n", n, ilsch->n);
    if (ilev) {
        data->schur_start = nn - ilsch->bsz[ilsch->n];//data->schur_start = n - ilsch->n;//ilsch->n is the order of C
        //data->schur_start = nn - ilsch->bsz[ilsch->n];//point index
    }
    else {
        is = self->is;
        data->schur_start = schur_start_p;
        //data->schur_start_p = schur_start_p;//point index
    }
    //  printf("data->schur_start = %d, data->schur_start_p = %d\n", data->schur_start, data->schur_start_p);
    PARMS_MEMCPY(data->ipar,   param->ipar, 18);
    PARMS_MEMCPY(data->pgfpar, param->pgfpar, 2);
    return 0;
}/*-----end-of-BARMS2----------------------------------------------------
   +--------------------------------------------------------------------*/

int parms_barms_vcsr_l(parms_Mat self, parms_FactParam param, void *mat,
                       parms_Operator *op)
{
    int nBlock, *nB = NULL, *perm = NULL, ierr, i, n, *iperm = NULL;
    double eps, blocksize, Bdensity;// tib,
    parms_vcsr csmat = NULL;;
    parms_bvcsr vbmat = NULL;
    parms_Map  is;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    printf("in parms_barms_vcsr_l\n");


    is = self->is;

    csmat = (parms_vcsr)mat;
    n = csmat->n;

    eps = param->eps;


//    outputcsmatpa(csmat,"csmatlocal",1);

    ierr = init_blocks( csmat, &nBlock, &nB, &perm, eps);//int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps)//parms_PCSetup(pc);

    if(ierr != 0) {
        //fprintf(stderr, "ierr = %d\n", ierr);
        fprintf(stderr, "*** in init_blocks ierr != 0 ***\n");
        //MPI_Finalize();
        exit(1);
    }
    //output_intvector("perm.coo",perm,0, n);getchar();
    if( dpermC( csmat, perm ) != 0 ) {
        fprintf( stderr, "*** dpermC error ***\n" );
        MPI_Finalize();
        exit(1);
    }

    param->n = nBlock;
    param->schur_start = param->n;
    //printf("(*op)->ref = %d \n",(*op)->ref);

    /* parms_Operator newOpt; */
    /* newOpt = *op; */
    /* newOpt->blockperm = perm; */// move permutation array to operator struct for solving phase
    //(*op)->blockperm = perm;// move permutation array to operator struct for solving phase

    //param->n = param->eps;//missing

    /*-------------------- convert to block matrix. */
    vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
    ierr = csrvbsrC( 1, nBlock, nB, csmat, vbmat );

    PARMS_NEWARRAY(iperm, n);
    /*-------------------- obtain reverse permutation array   */
    for (i=0; i<n; i++)
        iperm[perm[i]] = i;

    if( dpermC( csmat, iperm ) != 0 ) {
        fprintf( stderr, "*** dpermC error ***\n" );
        MPI_Finalize();
        exit(1);
    }

    //outputcsmatpa(csmat,"csmatipermed",1);
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1);  */
    PARMS_FREE(iperm);
    PARMS_FREE(nB);

    blocksize = (double)csmat->n / (double)nBlock;
    Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;

    //printf("is->schur_start = %d, blocksize = %f \n",is->schur_start, blocksize);

    is->schur_start = (int)(is->schur_start/blocksize);//just a approximate number of internal nodes
    //printf("is->schur_start = %d \n",is->schur_start);

    // if (myid == 0)
    printf("\n Bsize=%-7f, Bdensity=%-7f, myid = %d\n",blocksize, Bdensity, myid);
    //outputvbmatpa(vbmat,"vbmatlocal",1);

    parms_barms_vcsr(self, param, vbmat, op);

    (*op)->blockperm = perm;// move permutation array to operator struct for solving phase

    //outputvbmatpa(vbmat,"vbmatlocal",1);

    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */
    return 0;
}/*-----end-of-localgcBARMS2----------------------------------------------------
   +--------------------------------------------------------------------*/



int parms_barms_vcsr_old(parms_Mat self, parms_FactParam param, void *mat,
                         parms_Operator *op)
{
    /*---------------------------------------------------------------------
| MULTI-LEVEL BLOCK ILUT PRECONDITIONER.
| ealier  version:  June 23, 1999  BJS -- 
| version2: Dec. 07th, 2000, YS  [reorganized ]
| version 3 (latest) Aug. 2005.  [reorganized + includes ddpq]
+---------------------------------------------------------------------- 
| ON ENTRY:
| ========= 
| ( Amat ) = original matrix A stored in C-style Compressed Sparse
|            Row (CSR) format -- 
|            see LIB/heads.h for the formal definition of this format.
|
| ipar[0:17]  = integer array to store parameters for 
|       vbarms construction (arms2) 
|
|       ipar[0]:=nlev.  number of levels (reduction processes). 
|                       see also "on return" below. 
| 
|       ipar[1]:= level-reordering option to be used.  
|                 if ipar[1]==0 ARMS uses a block independent set ordering
|                  with a sort of reverse cutill Mc Kee ordering to build 
|                  the blocks. This yields a symmetric ordering. 
|                  in this case the reordering routine called is indsetC
|                 if ipar[1] == 1, then a nonsymmetric ordering is used.
|                  In this case, the B matrix is constructed to be as
|                  diagonally dominant as possible and as sparse as possble.
|                  in this case the reordering routine called is ddPQ.
|                 
|       ipar[2]:=bsize. for indset  Dimension of the blocks. 
|                  bsize is only a target block size. The size of 
|                  each block can vary and is >= bsize. 
|                  for ddPQ: this is only the smallest size of the 
|                  last level. vbarms will stop when either the number 
|                  of levels reaches nlev (ipar[0]) or the size of the
|                  next level (C block) is less than bsize.
|
|       ipar[3]:=iout   if (iout > 0) statistics on the run are 
|                       printed to FILE *ft
|
|	ipar[4]:= Krylov subspace dimension for last level 
|		    ipar[4] == 0 means only backward/forward solve
|		    is performed on last level.                  
|	ipar[5]:=  maximum # iterations on last level
|
|       ipar[6-9] NOT used [reserved for later use] - set to zero.
| 
| The following set method options for vbarms2. Their default values can
| all be set to zero if desired. 
|
|       ipar[10-13] == meth[0:3] = method flags for interlevel blocks
|       ipar[14-17] == meth[0:3] = method flags for last level block - 
|       with the following meaning in both cases:
|            meth[0] nonsummetric permutations of  1: yes. affects rperm
|                    USED FOR LAST SCHUR COMPLEMENT 
|            meth[1] permutations of columns 0:no 1: yes. So far this is
|                    USED ONLY FOR LAST BLOCK [ILUTP instead of ILUT]. 
|                    (so ipar[11] does no matter - enter zero). If 
|                    ipar[15] is one then ILUTP will be used instead 
|                    of ILUT. Permutation data stored in: perm2. 
|            meth[2] diag. row scaling. 0:no 1:yes. Data: D1
|            meth[3] diag. column scaling. 0:no 1:yes. Data: D2
|       all transformations related to parametres in meth[*] (permutation, 
|       scaling,..) are applied to the matrix before processing it 
| 
| ft       =  file for printing statistics on run
|
| droptol  = Threshold parameters for dropping elements in ILU 
|            factorization.
|            droptol[0:4] = related to the multilevel  block factorization
|            droptol[5:6] = related to ILU factorization of last block.
|            This flexibility is more than is really needed. one can use
|            a single parameter for all. it is preferable to use one value
|            for droptol[0:4] and another (smaller) for droptol[5:6]
|            droptol[0] = threshold for dropping  in L [B]. See piluNEW.c:
|            droptol[1] = threshold for dropping  in U [B].
|            droptol[2] = threshold for dropping  in L^{-1} F 
|            droptol[3] = threshold for dropping  in E U^{-1} 
|            droptol[4] = threshold for dropping  in Schur complement
|            droptol[5] = threshold for dropping  in L in last block
|              [see ilutpC.c]
|            droptol[6] = threshold for dropping  in U in last block
|              [see ilutpC.c]
|             This provides a rich selection - though in practice only 4
|             parameters are needed [which can be set to be the same 
              actually] -- indeed it makes sense to take
|             droptol[0] = droptol[1],  droptol[2] = droptol[3], 
|             and droptol[4] = droptol[5]
|
| lfil     = lfil[0:6] is an array containing the fill-in parameters.
|            similar explanations as above, namely:
|            lfil[0] = amount of fill-in kept  in L [B]. 
|            lfil[1] = amount of fill-in kept  in U [B].
|            lfil[2] = amount of fill-in kept  in E L\inv 
|            lfil[3] = amount of fill-in kept  in U \inv F
|            lfil[4] = amount of fill-in kept  in S    .
|            lfil[5] = amount of fill-in kept  in L_S  .
|            lfil[6] = amount of fill-in kept  in U_S 
|             
| tolind   = tolerance parameter used by the indset function. 
|            a row is not accepted into the independent set if 
|            the *relative* diagonal tolerance is below tolind.
|            see indset function for details. Good values are 
|            between 0.05 and 0.5 -- larger values tend to be better
|            for harder problems.
| 
| ON RETURN:
|=============
|
| (PreMat)  = vbarms data structure which consists of two parts:
|             levmat and ilsch. 
|
| ++(levmat)= permuted and sorted matrices for each level of the block 
|             factorization stored in PerMat4 struct. Specifically
|             each level in levmat contains the 4 matrices in:
|
|
|            |\         |       |
|            |  \   U   |       |
|            |    \     |   F   |
|            |  L   \   |       |
|            |        \ |       |
|            |----------+-------|
|            |          |       |
|            |    E     |       |
|            |          |       |
|            
|            plus a few other things. See LIB/heads.h for details.
|
| ++(ilsch) = IluSpar struct. If the block of the last level is:
|
|                        |  B    F |
|                  A_l = |         | 
|                        |  E    C |
|
|             then IluSpar contains the block C and an ILU
|             factorization (matrices L and U) for the last 
|             Schur complement [S ~ C - E inv(B) F ]
|             (modulo dropping) see LIB/heads.h for details.
|
| ipar[0]   = number of levels found (may differ from input value) 
|
+---------------------------------------------------------------------*/
    /* local matrix object */
    parms_Operator newOpt;
    parms_barms_data data;
    parms_bvcsr      Amat;
    parms_Map       is;
    /*-------------------- function  prototyping  done in LIB/protos.h    */
    /*-------------------- move above to protos.h */
    vbp4ptr levp=NULL, levc=NULL, levn=NULL, levmat=NULL;
    vbsptr schur=NULL, B=NULL, F=NULL, E=NULL, C=NULL;
    vbilutptr ilsch=NULL;
    csptr csmat2=NULL;
    /*-------------------- local variables  (initialized)   */
    double *dd1 = NULL, *dd2 = NULL;
    double *droptol, tolind;
    int nlev, bsize, iout, ierr = 0;
    int Storedd12 = 0; //for judging dd1 and dd2 should be stored or not
    int *ipar, *lfil;
    int methL[4], methS[4];
    /*--------------------  local variables  (not initialized)   */
    int nA, nB, nC, j, n, ilev, symperm, schur_start, nbnd, nn, schur_start_p;//i,
    FILE *ft;
    /*--------------------    work arrays:    */
    int *iwork, *uwork, *iworkn, *uworkn, *nBB, nset;
    /*   timer arrays:  */
    /*   double *symtime, *unstime, *factime, *tottime;*/
    /*---------------------------BEGIN ARMS-------------------------------*/
    /*   schur matrix starts being original A */

    /*-------------------- begin                                         */
    n = param->n;
    nbnd = schur_start = param->schur_start;
    is = self->is;
    if (schur_start == -1) {
        nbnd = is->schur_start;
        printf("nbnd = %d \n",nbnd);
    }



    // This statement is invalid not only here but also in arms2 ????
    if (!param->isalloc) {
        // printf("if (!param->isalloc)\n");

        parms_OperatorCreate(&newOpt);
        PARMS_MEMCPY(newOpt->ops, &parms_barms_sol_vptr, 1);
        PARMS_NEW(data);
        PARMS_NEW(data->levmat);
        levmat = data->levmat;
        //printf("data->levmat=%p\n", data->levmat);
        PARMS_NEW(data->ilus);
        ilsch =  data->ilus;
        newOpt->data = data;
        *op = newOpt;
    }
    else {
        //printf("in else\n");
        data = (*op)->data;
        //printf("nothing wrong\n");
    }
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */
    Amat = (parms_bvcsr)mat;

    /* compute the number of nonzero entries in block form mat */
    data->nnz_mat = nnzVBMat1(Amat);//memVBMat(Amat);//nnzVBMat1(Amat);

    /* int myid; */
    /* MPI_Comm_rank(MPI_COMM_WORLD, &myid); */
    /* printf("The number of non-zero entries is %d in local block form matrix, pid = %d \n",data->nnz_mat, myid); */
    /* printf("The number of stored entries is %d in local block form matrix, pid = %d \n",memVBMat(Amat), myid); */
    /* outputvbmatpa(Amat,"localvbmat_data",1); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */

    /* retrieve data from input objects */
    ipar = param->ipar;
    lfil = param->lfil;
    droptol = param->droptol;
    tolind  = param->tolind;

    nlev = ipar[0];
    //printf("nlev=%d\n", nlev);
//    printf("in barms_vcsr\n");
    fflush(stdout);
    bsize = ipar[2];
    iout = ipar[3];

    if (iout > 0 ) {
        ft = stdout;
    }
    else{
        ft = NULL;
    }

    ierr = 0;
    /*---------------------------------------------------------------------
| The matrix (a,ja,ia) plays role of Schur compl. from the 0th level.
+--------------------------------------------------------------------*/

    nC = nA = n = Amat->n;
    nn = Amat->bsz[n];

    schur_start_p = Amat->bsz[is->schur_start]; //new, for the point index of schur_start, needed in schur-solver

    //schur_start_p = Amat->bsz[is->schur_start]; //new, for the point index of schur_start, needed in schur-solver

    if (bsize >= n) bsize = n-1;
    levmat->n = n; levmat->nB = 0;

    //FIXME: Not sure why this is needed here and not in arms2
    //schur = (vbsptr) Malloc(sizeof(VBSparMat), "arms2:7" );
    //vbscpy(Amat, schur);
    schur = Amat;
    levc = levmat;
    /*--------------------------------------- */
    levc->prev = levc->next = levp = NULL;
    levc->n = 0;

    iwork = (int *) Malloc(n*sizeof(int), "barms2:2.1" );
    nBB = (int *) Malloc(n*sizeof(int), "barms2:2.2" );

    memcpy(methL, &ipar[10], 4*sizeof(int));
    memcpy(methS, &ipar[14], 4*sizeof(int));
    /*---------------------------------------------------------------------
    | The preconditioner construction is divided into two parts:
    |   1st part: construct and store multi-level L and U factors;
    |   2nd part: construct the ILUT factorization for the coarsest level
    +--------------------------------------------------------------------*/
    if ( (iout > 0)  && (nlev > 0) ) {
        fprintf(ft,"  \n");
        fprintf(ft,"Level   Total Unknowns    B-block   Coarse set\n");
        fprintf(ft,"=====   ==============    =======   ==========\n");
    }
    /*---------------------------------------------------------------------
    | main loop to construct multi-level LU preconditioner. Loop is on the
    | level ilev. nA is the dimension of matrix A_l at each level.
    +--------------------------------------------------------------------*/
    for (ilev = 0; ilev < nlev; ilev++) {
        /*-------------------- new nA is old nC -- from previous level */
        nA = nC;

        if ( nA <= bsize )  goto label1000;


        /*-------------------- allocate work space                        */
        iworkn = (int *) Malloc(schur->bsz[nA]*sizeof(int), "barms2:2.5" );
        symperm = 0;    /* 0nly needed in cleanP4 */
        if (ipar[1] == 1) {
            uworkn = (int *) Malloc(schur->bsz[nA]*sizeof(int), "barms2:2.6" );
            uwork = (int *) Malloc(n*sizeof(int), "barms2:2.6" );
        }
        else{
            symperm = 1;
            uwork = iwork;
            uworkn = iworkn;
        }
        /*-------------------- SCALING*/
        dd1 = NULL;
        dd2 = NULL;
        if (methL[2]) {
            dd1 = (double *) Malloc(schur->bsz[nA]*sizeof(double), "barms2:3" );
            j=vbfroscalC(schur, dd1,1);
            if (j) printf("ERROR in vbfroscalC -  row %d  is a zero row\n",j);
        }

        if (methL[3]) {
            dd2 = (double *) Malloc(schur->bsz[nA]*sizeof(double), "barms2:4" );
            j=vbfcoscalC(schur, dd2,1);
            if (j) printf("ERROR in vbfcoscalC - column %d is a zero column\n",j);
        }
        /*--------------------independent-sets-permutation-------------------
      |  do reordering -- The matrix and its transpose are used.
      +--------------------------------------------------------------------*/
        /* if (SHIFTTOL > 0.0) shiftsD(schur,SHIFTTOL);    */
        csmat2 = (csptr) malloc(sizeof(SparMat));
        ierr = vbsrc2csr(schur, csmat2);
        if (ipar[1] == 1)
            PQperm(csmat2, uwork, iwork, &nB, tolind, nbnd) ;
        else
            //  indsetC(csmat2, bsize, iwork, &nB, tolind, nbnd) ;
            indsetC2(csmat2, bsize, iwork, &nB, tolind, nbnd, nBB, &nset) ;
        cleanCS(csmat2);
        /*---------------------------------------------------------------------
      | nB is the total number of nodes in the independent set.
      | nC : nA - nB = the size of the reduced system.
      +--------------------------------------------------------------------*/
        nC = nA - nB;
        nbnd -= nB;
        /*   if the size of B or C is zero , exit the main loop  */
        /*   printf ("  nB %d nC %d \n",nB, nC); */
        //printf("nB=%d\n", nB);
        //printf("nC=%d\n", nC);
        if ( nB == 0 || nC == 0 )  {
            Storedd12 = 1;
            free(iworkn);iworkn = NULL;
            goto label1000;
        }
        /*int myid;
MPI_Comm_rank(MPI_COMM_WORLD, &myid);
printf("myid=%d\n", ilev);
  printf("ilev=%d\n", ilev);
  printf("nA=%d\n", nA);
  printf("bsize=%d\n", bsize);*/
        /*---------------------------------------------------------------------
      | The matrix for the current level is in (schur).
      | The permutations arrays are in iwork and uwork (row).
      | The routines rpermC, cpermC permute the matrix in place.
      *-----------------------------------------------------------------------*/
        /*   DEBUG : SHOULD THIS GO BEFORE GOTO LABEL1000 ?? */
        vbrpermC(schur,uwork);
        vbcpermC(schur,iwork);
        ierr = permuten(iwork, iworkn, schur->bsz, nA);
        /*   prtC(schur, ilev) ;   print matrix - debugging */
        /*-----------------------------------------------------------------------
      | If this is the first level, the permuted matrix is stored in
      | (levc) = (levmat).  Otherwise, the next level is created in (levc).
      +--------------------------------------------------------------------*/
        if (ilev > 0) {
            /*-   delete C matrix of any level except last one (no longer needed) */
            cleanVBMat(C);
            levn = (vbp4ptr) Malloc(sizeof(VBPer4Mat), "barms2:6" );
            /* levc->prev = levp; */
            levc->next = levn;
            levp = levc;
            levc = levn;
            levc->prev = levp;
        }
        /*-------------------- p4ptr struct for current schur complement */
        B = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:7" );
        E = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:8" );
        F = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:9" );
        C = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:10" );
        vbSplit4copy(schur, nB, nC, B, F, E, C);
        ierr = setupVBP4_vbarmsold(levc, nB, nC, F, E, schur->bsz);
//        printf("levc->lu value is %p, in barms.\n", levc->lu);//%f %p %s %c

        /*--------------------     copy a few pointers       ---- */
        levc->perm  = iworkn;
        levc->rperm = uworkn;
        levc->symperm = symperm;
        levc->D1=dd1;
        levc->D2=dd2;
        /*---------------------------------------------------------------------
      | a copy of the matrix (schur) has been permuted. Now perform the
      | block factorization:
      |
      | | B   F |       | L       0 |     | U  L^-1 F |
      | |       |   =   |           |  X  |           | = L x U
      | | E   C |       | E U^-1  I |     | 0    A1   |
      |
      | The factors E U^-1 and L^-1 F are discarded after the factorization.
      |
      +--------------------------------------------------------------------*/
        if (iout > 0)
            fprintf(ft,"%3d %13d %13d %10d\n", ilev+1,schur->bsz[nA],schur->bsz[nB],schur->bsz[nA]-schur->bsz[nB]);
        /*---------------------------------------------------------------------
      | PILUT constructs one level of the block ILU fact.  The permuted matrix
      | is in (levc).  The L and U factors will be stored in the p4mat struct.
      | destroy current Schur  complement - no longer needed  - and set-up new
      | one for next level...
      +--------------------------------------------------------------------*/
        cleanVBMat(schur);
        schur = (vbsptr) Malloc(sizeof(VBSparMat), "barms2:11" );
        /*----------------------------------------------------------------------
      | calling PILU to construct this level block factorization
      | ! core dump in extreme case of empty matrices.
      +----------------------------------------------------------------------*/
        //FIXME!!!!!
        //~ if (ilutp == 1)
        //~ ierr = computschurpartitionp(levc, B, C, droptol, lfil, ipar, schur, nBB, nset);
        //~ else
        //~ printf("levc=%p\n", levc);
        ierr = computschurpartition(levc, B, C, droptol, lfil, schur, nBB, nset);
        //~ setupVBMat(schur, nC, NULL);
        //~ ierr = vbiluNEW(levc, B, C, droptol, lfil, schur);
        //~ outputvbmat(schur, "newschur", 1);
        //~ outputvbmat(C, "cvbarms", 1);
        //ierr = computschurpartition(levc, B, C, droptol, lfil, schur, nBB, nset);
        //~ outputvbmat(schur, "oldschur", 1);
        /* prtC(levc->L, ilev) ; */
        if (ierr && ft) {
            fprintf(ft," ERROR IN  PILU  -- IERR = %d\n", ierr);
            return(1);
        }
        cleanVBMat(B);
    }
    /*---------------------------------------------------------------------
    |   done with the reduction. Record the number of levels in ipar[0]
    |**********************************************************************
    +--------------------------------------------------------------------*/
label1000:
    /* printf (" nnz_Schur %d \n",cs_nnz (schur)); */
    if (iwork) free(iwork);
    if (nBB) free(nBB);
    levc->next = NULL;
    ipar[0] = ilev;
    data->nlev = ilev;
    data->n = nn;
    nC = schur->n;
    setupVBILUT(ilsch, nC, schur->bsz);

    if ( Storedd12 )  {
        ilsch->dd1 = dd1;
        ilsch->dd2 = dd2;
    }
    /*--------------------------------------------------------------------*/
    /* define C-matrix (member of ilsch) to be last C matrix */
    if (ilev > 0) ilsch->C=C;

    /*-------------------- for ilut fact of schur matrix */
    /*  SCALING  */

    ilsch->D1 = NULL;
    if (methS[2] || methS[1] == 0) {
        //printf("in scale\n");
        ilsch->D1 = (double *) Malloc(schur->bsz[nC]*sizeof(double), "barms2:iluschD1" );
        j=vbfroscalC(schur, ilsch->D1, 1);
        if (j) printf("ERROR in vbfroscalC - row %d is a zero row\n",j);
    }

    ilsch->D2  = NULL;
    if (methS[3]) {
        ilsch->D2 = (double *) Malloc(schur->bsz[nC]*sizeof(double), "barms2:iluschD1" );
        j =vbfcoscalC(schur, ilsch->D2, 1);
        if (j) printf("ERROR in vbfcoscalC - column %d is a zero column\n",j);
    }

    /*---------------------------------------------------------------------
    |     get ILUT factorization for the last reduced system.
    +--------------------------------------------------------------------*/
    //~ uwork = NULL;
    iwork = NULL;

    ilsch->rperm = iwork; //uwork;
    ilsch->perm  = iwork;

    /*   printf("  lf : %d  %d  %d  %d  %d  %d  %d  \n",lfil[0],
       lfil[1], lfil[2], lfil[3], lfil[4], lfil[5], lfil[6]) ; */

    ierr = vbilutD(schur, droptol, lfil, ilsch);
    /*---------- OPTIMIZATRION: NEED TO COMPOUND THE TWO
    RIGHT PERMUTATIONS -- CHANGES HERE AND IN
    USCHUR SOLVE ==  compound permutations */
    if (ierr && ft) {
        fprintf(ft," ERROR IN  ILUT -- IERR = %d\n", ierr);
        return(1);
    }
    /* Last Schur complement no longer needed */
    cleanVBMat(schur);
    data->nnz_prec = nnz_vbarms(data, stderr);
    printf("bsz %d %d\n", ilsch->bsz[ilsch->n], ilsch->lu->bsz[ilsch->lu->n]);
    //~ data->ind = n - ilsch->n;
    data->ind = nn - ilsch->bsz[ilsch->n];
    //  printf("n=%d, ilschn=%d\n", n, ilsch->n);
    if (ilev) {
        data->schur_start = nn - ilsch->bsz[ilsch->n];//data->schur_start = n - ilsch->n;//ilsch->n is the order of C
        //data->schur_start = nn - ilsch->bsz[ilsch->n];//point index
    }
    else {
        is = self->is;
        data->schur_start = schur_start_p;
        //data->schur_start_p = schur_start_p;//point index
    }
    //  printf("data->schur_start = %d, data->schur_start_p = %d\n", data->schur_start, data->schur_start_p);
    PARMS_MEMCPY(data->ipar,   param->ipar, 18);
    PARMS_MEMCPY(data->pgfpar, param->pgfpar, 2);
    return 0;
}/*-----end-of-BARMS2----------------------------------------------------
   +--------------------------------------------------------------------*/

int parms_barms_vcsr_l_old(parms_Mat self, parms_FactParam param, void *mat,
                           parms_Operator *op)
{
    int nBlock, *nB = NULL, *perm = NULL, ierr, i, n, *iperm = NULL;
    double eps, blocksize, Bdensity;// tib,
    parms_vcsr csmat = NULL;;
    parms_bvcsr vbmat = NULL;
    parms_Map  is;

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    printf("in parms_barms_vcsr_l\n");


    is = self->is;

    csmat = (parms_vcsr)mat;
    n = csmat->n;

    eps = param->eps;


//    outputcsmatpa(csmat,"csmatlocal",1);

    ierr = init_blocks( csmat, &nBlock, &nB, &perm, eps);//int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps)//parms_PCSetup(pc);

    if(ierr != 0) {
        //fprintf(stderr, "ierr = %d\n", ierr);
        fprintf(stderr, "*** in init_blocks ierr != 0 ***\n");
        //MPI_Finalize();
        exit(1);
    }
    //output_intvector("perm.coo",perm,0, n);getchar();
    if( dpermC( csmat, perm ) != 0 ) {
        fprintf( stderr, "*** dpermC error ***\n" );
        MPI_Finalize();
        exit(1);
    }

    param->n = nBlock;
    param->schur_start = param->n;
    //printf("(*op)->ref = %d \n",(*op)->ref);

    /* parms_Operator newOpt; */
    /* newOpt = *op; */
    /* newOpt->blockperm = perm; */// move permutation array to operator struct for solving phase
    //(*op)->blockperm = perm;// move permutation array to operator struct for solving phase

    //param->n = param->eps;//missing

    /*-------------------- convert to block matrix. */
    vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
    ierr = csrvbsrC( 1, nBlock, nB, csmat, vbmat );

    PARMS_NEWARRAY(iperm, n);
    /*-------------------- obtain reverse permutation array   */
    for (i=0; i<n; i++)
        iperm[perm[i]] = i;

    if( dpermC( csmat, iperm ) != 0 ) {
        fprintf( stderr, "*** dpermC error ***\n" );
        MPI_Finalize();
        exit(1);
    }

    //outputcsmatpa(csmat,"csmatipermed",1);
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1);  */
    PARMS_FREE(iperm);
    PARMS_FREE(nB);

    blocksize = (double)csmat->n / (double)nBlock;
    Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;

    //printf("is->schur_start = %d, blocksize = %f \n",is->schur_start, blocksize);

    is->schur_start = (int)(is->schur_start/blocksize);//just a approximate number of internal nodes
    //printf("is->schur_start = %d \n",is->schur_start);

    // if (myid == 0)
    printf("\n Bsize=%-7f, Bdensity=%-7f, myid = %d\n",blocksize, Bdensity, myid);
    //outputvbmatpa(vbmat,"vbmatlocal",1);

    parms_barms_vcsr_old(self, param, vbmat, op);

    (*op)->blockperm = perm;// move permutation array to operator struct for solving phase

    //outputvbmatpa(vbmat,"vbmatlocal",1);

    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */
    return 0;
}/*-----end-of-localgcBARMS2----------------------------------------------------
   +--------------------------------------------------------------------*/
