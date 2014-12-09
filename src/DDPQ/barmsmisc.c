#include "protos.h"

int setupVBILU(vbiluptr lu, int n, int *bsz);
int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol,
             BData *w, FILE *fp );

double vbnorm2( int sz, FLOAT *a )
{
/*-------------------- return average norm among a[sz] */
    int tmp = 1;

    return GNRM2( sz, a, tmp ) / (double)sz;
}

int cleanVBMat( vbsptr vbmat )
{
/*----------------------------------------------------------------------
| Free up memory allocated for VBSpaFmt structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat )  =  Pointer to a VBSpaFmt struct.
|--------------------------------------------------------------------*/
    int i, j; 
    if( vbmat == NULL ) return 0;
    if( vbmat->n < 1 ) return 0;
    
    for( i = 0; i < vbmat->n; i++ ) {
        if( vbmat->nzcount[i] > 0 ) {
            free( vbmat->ja[i] );
            if( vbmat->ba && vbmat->ba[i] ) {
                for( j = 0; j < vbmat->nzcount[i]; j++ ) {
                    free( vbmat->ba[i][j] );
                }
                free( vbmat->ba[i] );
            }
        }
        if( vbmat->D && vbmat->D[i] ) free( vbmat->D[i] );
    }
    if( vbmat->D ) free( vbmat->D );
    free( vbmat->ja );
    if( vbmat->ba ) free( vbmat->ba );
    free( vbmat->nzcount );
    if( vbmat->bsz ) free( vbmat->bsz );
    free( vbmat );
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBMat
|--------------------------------------------------------------------*/
int cleanVBMat1( vbsptr vbmat )
{
/*----------------------------------------------------------------------
| Free up memory allocated for VBSpaFmt structs, for exculsively rectangular matrix.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat )  =  Pointer to a VBSpaFmt struct.
|--------------------------------------------------------------------*/
    int i, j; 
    if( vbmat == NULL ) return 0;
    if( vbmat->n < 1 ) return 0;
    
    for( i = 0; i < vbmat->n; i++ ) {
        if( vbmat->nzcount[i] > 0 ) {
            free( vbmat->ja[i] );
            if( vbmat->ba && vbmat->ba[i] ) {
                for( j = 0; j < vbmat->nzcount[i]; j++ ) {
                    free( vbmat->ba[i][j] );
                }
                free( vbmat->ba[i] );
            }
        }
        if( vbmat->D && vbmat->D[i] ) free( vbmat->D[i] );
    }
    if( vbmat->D ) free( vbmat->D );
    free( vbmat->ja );
    if( vbmat->ba ) free( vbmat->ba );
    free( vbmat->nzcount );
    if( vbmat->bsz ) free( vbmat->bsz );
    if( vbmat->bszc ) free( vbmat->bszc );
    free( vbmat );
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBMat1
|--------------------------------------------------------------------*/


int cleanVBILU( vbiluptr lu ){
/*----------------------------------------------------------------------
| Free up memory allocated for VBILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
|   ( lu )  =  Pointer to a VBILUSpar struct.
|--------------------------------------------------------------------*/
    int n = lu->n, i;
    if( NULL == lu ) return 0;
    if( lu->D ) {
        for( i = 0; i < n; i++ ) {
            if( lu->D[i] ) free( lu->D[i] );
        }
        free( lu->D );
    }
    if( lu->bsz ) free( lu->bsz );
    cleanVBMat( lu->L );
    cleanVBMat( lu->U );
    if( lu->work ) free( lu->work );
    if( lu->bf ) free( lu->bf );
    free( lu );
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBILU
|--------------------------------------------------------------------*/

int cleanIILUT( iilutptr vbbmat ){
/*----------------------------------------------------------------------
| Free up memory allocated for VBILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|   ( vbbmat )  =  Pointer to a VBILUSpar struct.
|--------------------------------------------------------------------*/

//typedef struct IILUTfac {//for coarse point ilutp struct
    //int n;        /* the block dimension of the independent sets*/
    //int *bsz;     /* the row/col of the first element of each coarse block*/                        
    //ilutptr *llu;   /* pointer-to-pointer to store coarse blocks */
//} IIluSpar, *iilutptr;

/*  if (amat->lu) {
    cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    amat->lu = NULL;*/

    int n = vbbmat->n, i;
    if( NULL == vbbmat ) return 0;
    if( vbbmat->llu ) {
        for( i = 0; i < n; i++ ) {
            if( vbbmat->llu[i] ) {cleanILUT(vbbmat->llu[i], 0);free(vbbmat->llu[i]);vbbmat->llu[i]=NULL; }//int cleanILUT(ilutptr amat, int indic)
        }
        free( vbbmat->llu );
    }
    if( vbbmat->bsz ) free( vbbmat->bsz );
    free( vbbmat );
    return 0;
}



int vbsrc2csr(vbsptr vbmat, csptr csmat)
{
/*
 * This  subroutine converts the vbsptr matrix to csptr format, in which 
 * a(i,j)=norm|A(i,j)|, 
 * Here we choose Euclidean norm.  
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 *
 * vbmat = Various Block Sparse Row format Matrix, vbsptr
 *
 * on return:
 *-----------
 * 
 *  csmat = Sparse Row format Matrix, csptr
 *
 *  ierr  = integer, error code. 
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
  int i, j, nzcount, col, inc = 1, dimR, dimC, blocksz; 
  int n = vbmat->n, *ja, *bsz = vbmat->bsz;
  FLOAT *w;
  BData *ba;
  
  setupCS(csmat,n,1);
  w =(FLOAT*)malloc(n*sizeof(FLOAT));
  for( i = 0; i < n; i++ ) 
  {
      dimR = B_DIM(bsz,i);//get the row dimension of the block
      nzcount=csmat->nnzrow[i]=vbmat->nzcount[i];
      ja=vbmat->ja[i];
      csmat->pj[i] = (int *) Malloc(nzcount*sizeof(int), "vbsrc2csr:1" );
      memcpy(csmat->pj[i],ja,nzcount*sizeof(int));
      ba=vbmat->ba[i];
      for( j = 0; j < nzcount; j++ ) 
      {
          col = ja[j];
          dimC = B_DIM(bsz,col);//obtain the column dimension of the block
          blocksz = dimR*dimC;//caculate the blocksize
          w[j] = GNRM2(blocksz,ba[j],inc);/*caculate the Euclidean norm of the//FLOAT vbnorm2( int sz, FLOAT *a )//w[j]=vbnorm2( blocksz,ba[j]);
          //block and ready to be copied to a(i,j)*/
      }
      csmat->pa[i] = (FLOAT*)malloc(nzcount*sizeof(FLOAT));
      memcpy(csmat->pa[i],w,nzcount*sizeof(FLOAT));
  }
  free(w);
  return 0;
}

int vblusolC( FLOAT *y, FLOAT *x, vbiluptr lu)  
{
/*----------------------------------------------------------------------
 *    performs a forward followed by a backward block solve
 *    for LU matrix as produced by VBILUT
 *    y  = right-hand-side 
 *    x  = solution on return 
 *    lu = LU matrix as produced by VBILUT
 *
 *    note: lu->bf is used to store vector
 *--------------------------------------------------------------------*/
    int n = lu->n, *bsz = lu->bsz, i, j, bi, icol, dim, sz;
    int nzcount, nBs, nID, *ja, inc = 1, OPT;
    FLOAT *data, alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0;
    vbsptr L, U;
    BData *D, *ba;

    L = lu->L;
    U = lu->U;
    D = lu->D;
    OPT = lu->DiagOpt;
    /* Block L solve */
    for( i = 0; i < n; i++ ) {
        dim = B_DIM(bsz,i);
        nBs = bsz[i];
        for( j = 0; j < dim; j++ ) {
            nID = nBs + j;
            x[nID] = y[nID];
        }

        nzcount = L->nzcount[i];
        ja = L->ja[i];
        ba = L->ba[i];
        for( j = 0; j < nzcount; j++ ) {
            icol = ja[j];
            sz = B_DIM(bsz,icol);
            data = ba[j];
            GGEMV( "n",  dim,  sz,  alpha, data, dim, x+bsz[icol],
           inc, beta, x+nBs, inc );//GGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)//y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y
        }
    }
    /* Block -- U solve */
    for( i = n-1; i >= 0; i-- ) {
        dim = B_DIM(bsz,i);
        nzcount = U->nzcount[i];
        nBs = bsz[i];
        ja = U->ja[i];
        ba = U->ba[i];
        for( j = 0; j < nzcount; j++ ) {
            icol = ja[j];
            sz = B_DIM(bsz,icol);
            data = ba[j];
            GGEMV( "n", dim, sz, alpha, data, dim, x+bsz[icol], inc,
		   beta, x+nBs, inc ); 
        }
        data = D[i];
	if (OPT == 1) 
	  luinv( dim, data, x+nBs, lu->bf );
	else
      GGEMV( "n", dim, dim, alpha2, data, dim, x+nBs, inc, beta2,
         lu->bf, inc ); //GGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)//y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y
	
        for( bi = 0; bi < dim; bi++ ) {
            x[nBs+bi] = lu->bf[bi];
        }
    }

    return 0;
}

int vbrpermC(vbsptr vbmat, int *perm)
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the block rows of a matrix in VBSpaFmt format. 
| rperm  computes B = P A  where P is a permutation matrix.  
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination row number of row number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (amat) = a matrix stored in SpaFmt format.
|
|
| on return:
| ----------
| (amat) = P A stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
   int **addj, *nnz, i, size=vbmat->n, *nB, *addb;
   BData **addm;
   nB = (int *) Malloc( size*sizeof(int), "vbrpermC" );
   addj = (int **)Malloc( size*sizeof(int *), "vbrpermC" );
   addm = (BData **) Malloc( size*sizeof(BData *), "vbrpermC" );
   nnz = (int *) Malloc( size*sizeof(int), "vbrpermC" );
   addb = (int *) Malloc( size*sizeof(int), "vbrpermC" );
   for (i=0; i<size; i++) {
       nB[i] = B_DIM(vbmat->bsz,i);
       addb[perm[i]] = nB[i];
       addj[perm[i]] = vbmat->ja[i];
       addm[perm[i]] = vbmat->ba[i];
       nnz[perm[i]] = vbmat->nzcount[i];
   }
   for (i=0; i<size; i++) {
       if(i!=0){
       vbmat->bsz[i] = vbmat->bsz[i-1] + addb[i-1];
       }
       vbmat->ja[i] = addj[i];
       vbmat->ba[i] = addm[i];
       vbmat->nzcount[i] = nnz[i];
   }
   free(addj);
   free(addm);
   free(nnz);
   free(addb);
   free(nB);
   return 0;
}

int vbcpermC(vbsptr vbmat, int *perm) 
{
/*----------------------------------------------------------------------
|
| This subroutine permutes the columns of a matrix in SpaFmt format.
| cperm computes B = A P, where P is a permutation matrix.
| that maps column j into column perm(j), i.e., on return 
| The permutation P is defined through the array perm: for each j, 
| perm[j] represents the destination column number of column number j. 
|
|-----------------------------------------------------------------------
| on entry:
|----------
| (mat) = a matrix stored in SpaFmt format.
|
|
| on return:
| ----------
| (mat) = A P stored in SpaFmt format.
|
| integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|---------------------------------------------------------------------*/
   int i, j, *newj, size=vbmat->n, *aja;

	//printf("\ncperm row no= %4d \n",size); 

   newj = (int *) Malloc( size*sizeof(int), "vbcpermC" );
   for (i=0; i<size; i++) {
      aja = vbmat->ja[i];
      for (j=0; j<vbmat->nzcount[i]; j++)
	 newj[j] = perm[aja[j]];
  
      for (j=0; j<vbmat->nzcount[i]; j++)
	 aja[j] = newj[j];
      vbmat->ja[i] = aja;
   }
   free(newj);

   return 0;
}
int vbSplit4copy(vbsptr vbmat, int bsize, int csize, vbsptr B, vbsptr F,
	     vbsptr E, vbsptr C)
{
/*---------------------------------------------------------------------
| Convert permuted vbspaFmt struct to VBPerMat4 struct 
|                - matrix already permuted
|----------------------------------------------------------------------
| on entry:
|========== 
| ( vbmat )  =  Matrix stored in vbspaFmt format.
|              Internal pointers (and associated memory) destroyed before
|              return.
|
| On return:
|===========
|
| B, E, F, C = 4 blocks in 
| 
|          | B   F |      
|   vbmat= |       | 
|          | E   C | 
| 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int j, j1, numr, numl, ind, newj, rowz, *rowj, *new1j, *new2j, *nB, *nC, dimR, dimC;
   int *bsz = vbmat->bsz;
   BData *rowm, *new1m, *new2m;
/*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/
   nB = (int *) Malloc(bsize*sizeof(int), "vbSplit4copy:0" );
   for (j=0; j<bsize; j++){
      nB[j]=B_DIM(bsz,j);
   }

   nC = (int *) Malloc(csize*sizeof(int), "vbSplit4copy:0" );   
   for (j=0; j<csize; j++){
      nC[j]=B_DIM(bsz,j+bsize);
   }
   //output_intvector("nB.coo" ,nB,0, bsize);
   //output_intvector("nC.coo" ,nC,0, csize);		
   //printf("\n bsize=%d,csize=%d",bsize,csize); getchar();   
   if (setupVBMat(B,bsize,nB)) goto label111; 
   if (setupVBMat1(F,bsize,csize,nB,nC)) goto label111;
   if (setupVBMat1(E,csize,bsize,nC,nB)) goto label111; 
 //int setupVBMat1( vbsptr vbmat, int n,int nc, int *nBR, int *nBC )
 
   //output_intvector("nB.coo" ,nB,0, bsize);
   //output_intvector("nC.coo" ,nC,0, csize);		
   //output_intvector("FBSZ.coo" ,F->bsz,0, bsize+1);
   //output_intvector("FBSZC.coo" ,F->bszc,0, csize+1);
   //output_intvector("EBSZ.coo" ,E->bsz,0, csize+1);
   //output_intvector("EBSZC.coo" ,E->bszc,0, bsize+1);  
   //printf("\n Fbsize=%d,Fcsize=%d",F->n,F->nc); getchar();
   //printf("\n Ebsize=%d,Ecsize=%d",E->n,E->nc); getchar();
   if (setupVBMat(C,csize,nC)) goto label111;
	//~ output_intvector("nC.coo" ,nC,0, csize);	
   new1j = (int *) Malloc(bsize*sizeof(int), "vbSplit4copy:1" );
   new2j = (int *) Malloc(csize*sizeof(int), "vbSplit4copy:2" );
   new1m = (BData *) Malloc(bsize*sizeof(BData), "vbSplit4copy:3" );
   new2m = (BData *) Malloc(csize*sizeof(BData), "vbSplit4copy:4" );
/*    B and F blocks */ 
   for (j=0; j<bsize; j++) {
      numl = numr = 0;
      dimR = B_DIM(bsz,j);
      rowz = vbmat->nzcount[j];
      rowj = vbmat->ja[j];
      rowm = vbmat->ba[j];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      B->nzcount[j] = numl;
      F->nzcount[j] = numr;
      if (numl>0) {
	 B->ja[j] = (int *) Malloc(numl*sizeof(int), "vbSplit4copy:5" );
	 B->ba[j] = (BData *) Malloc(numl*sizeof(BData), "vbSplit4copy:6" );
      }
      if (numr>0) {
	 F->ja[j] = (int *) Malloc(numr*sizeof(int), "vbSplit4copy:7" );
	 F->ba[j] = (BData *) Malloc(numr*sizeof(BData), "vbSplit4copy:8" );
      }
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	 newj = rowj[j1];
	 dimC = B_DIM(bsz,newj);
	 if (newj<bsize) {
	    new1j[numl] = newj;
	    new1m[numl] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
	    copyBData( dimC, dimR, new1m[numl], rowm[j1], 0 );//new1m[numl] = rowm[j1];//??
	    numl++;
	 }
	 else {
	    new2j[numr] = newj - bsize;
	    new2m[numr] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
	    copyBData( dimC, dimR, new2m[numr], rowm[j1], 0 );//new2m[numr] = rowm[j1];//??
	    numr++;
	 }
      }
      memcpy(B->ja[j], new1j, numl*sizeof(int));
      memcpy(B->ba[j], new1m, numl*sizeof(BData));//??
      memcpy(F->ja[j], new2j, numr*sizeof(int));
      memcpy(F->ba[j], new2m, numr*sizeof(BData));
   }
/*    E and C blocks */
   for (j=0; j<csize; j++) {
      numl = numr = 0;
      ind = bsize + j;
      dimR = B_DIM(bsz,ind);
      rowz = vbmat->nzcount[ind];
      rowj = vbmat->ja[ind];
      rowm = vbmat->ba[ind];
      for (j1=0; j1<rowz; j1++) {
	 if (rowj[j1]<bsize) numl++;
	 else numr++;
      }
      E->nzcount[j] = numl;
      C->nzcount[j] = numr;
      if (numl>0) {
	E->ja[j] = (int *) Malloc(numl*sizeof(int), "vbSplit4copy:9" );
	E->ba[j] = (BData *) Malloc(numl*sizeof(BData), "vbSplit4copy:10" );
      }	
      if (numr>0) {
	C->ja[j] = (int *) Malloc(numr*sizeof(int), "vbSplit4copy:11" );
	C->ba[j] = (BData *) Malloc(numr*sizeof(BData), "vbSplit4copy:12" );
      }		
      numl = numr = 0;
      for (j1=0; j1<rowz; j1++) {
	newj = rowj[j1];
	dimC = B_DIM(bsz,newj);
	if (newj<bsize) {
	  new1j[numl] = newj;
	  new1m[numl] = (FLOAT *) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
	  copyBData( dimC, dimR, new1m[numl], rowm[j1], 0 );	  //new1m[numl] = rowm[j1];
	  numl++;
	}
	else {
	  new2j[numr] = newj - bsize;
	  new2m[numr] = (FLOAT *) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
	  copyBData( dimC, dimR, new2m[numr], rowm[j1], 0 );	  //new2m[numr] = rowm[j1];
	  numr++;
	}
      }
      memcpy(E->ja[j], new1j, numl*sizeof(int));
      memcpy(E->ba[j], new1m, numl*sizeof(BData));
      memcpy(C->ja[j], new2j, numr*sizeof(int));
      memcpy(C->ba[j], new2m, numr*sizeof(BData));
   }

   if (new1j) free(new1j);
   if (new2j) free(new2j);
   if (new1m) free(new1m);
   if (new2m) free(new2m);
   if (nB) free(nB);
   if (nC) free(nC);	
   return 0;
label111:
   return 1;
}

int setupVBP4 (vbp4ptr vbmat, int Bn, int Cn,  vbsptr F,  vbsptr E , int *bsz) // not yet
{
/*----------------------------------------------------------------------
| initialize VBPerMat4 struct given the F, E, blocks.  
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat )  =  Pointer to a PerMat4 struct.
|     Bn    =  block dimension of B block
|     Cn    =  block dimension of C block
|     F, E  = the two blocks to be assigned to struct - without the
|
| On return:
|===========
|
|  vbmat->n
|	->bsz
|	->lu                for each block: vbmat->M->n
|                                             ->nzcount
|      ->E                                       ->ja
|      ->F                                       ->ma
|      ->perm
|      ->rperm       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int n,i,nn;
   //int max_block_size = sizeof(FLOAT)*MAX_BLOCK_SIZE*MAX_BLOCK_SIZE;
   /* size n */
   n = vbmat->n = Bn + Cn;
   nn = bsz[n];
   vbmat->bsz=(int *) Malloc((n+1)*sizeof(int), "setupVBP4:2" );
   for( i = 0; i <= n; i++ ) vbmat->bsz[i] = bsz[i];//??
//printf("\n vbmat->bsz=%d,max=%d",vbmat->bsz,max_block_size); getchar();
   vbmat->nB = Bn;
   //nnB=(int *) Malloc(Bn*sizeof(int), "setupVBP4:2" );
   //for(i=0;i<Bn;i++){
   //nnB[i]=B_DIM(bsz,i);
//printf("\n i=%d,nnb[i]=%d",i,nnB[i]); getchar();
/* vbmat->perm = (int *) Malloc(n*sizeof(int), "setupP4:1" ); */
/*   assign space for wk -- note that this is only done at 1st level
     at other levels, copy pointer of wk from previous level */
   if (vbmat->prev == NULL)  /* wk has 2 * n entries now */
     vbmat->wk   = (FLOAT *) Malloc(2*nn*sizeof(FLOAT), "setupVBP4:2" );
   else 
     vbmat->wk = (vbmat->prev)->wk; 

/*-------------------- L and U */ 
   //vbmat->L = (vbsptr) Malloc(sizeof(VBSparMat), "setupVBP4:3" );
   //if (setupVBMat(vbmat->L, Bn,nnB)) return 1;//int setupVBMat( vbsptr vbmat, int n, int *nB )
   /*    fprintf(stdout,"  -- BN %d   Cn   %d \n", Bn,Cn);  */
   //vbmat->U = (vbsptr) Malloc(sizeof(VBSparMat), "setupVBP4:4" );
   //if (setupVBMat(vbmat->U, Bn,nnB)) return 1;

   // new sven
   vbmat->lu = (vbiluptr) Malloc(sizeof(VBILUSpar), "setupVBP4:3" );
   setupVBILU(vbmat->lu, Bn, bsz);

   //~ vbmat->lu = NULL;
   vbmat->F = F; 
   vbmat->E = E; 
   vbmat->plu = NULL;
   //vbmat->work = (int *)Malloc( sizeof(int) * n, "setupVBILU" );
   //vbmat->bf = (BData)Malloc( max_block_size, "setupVBILU" );
   return 0;
}
/*---------------------------------------------------------------------
|     end of setupVBP4 
|--------------------------------------------------------------------*/
int setupVBILUT(vbilutptr vbmat, int len, int *bsz)
{
/*----------------------------------------------------------------------
| Allocate pointers for ILUTfac structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbmat )  =  Pointer to a ILUTfac struct.
|     len   =  size of L U  blocks
|
| On return:
|===========
|
|  vbmat->L                for each block: vbmat->M->n
|      ->U                                       ->nzcount
|                                                ->ja
|                                                ->ma
|      ->rperm       (if meth[0] > 0)
|      ->perm2       (if meth[1] > 0)
|      ->D1          (if meth[2] > 0)
|      ->D2          (if meth[3] > 0)
|
|  Permutation arrays are initialized to the identity.
|  Scaling arrays are initialized to 1.0.
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
  int n,i,nn;
  n = vbmat->n = len;
  nn = bsz[n];
  vbmat->bsz=(int *) Malloc((n+1)*sizeof(int), "setupVBILUT:2" );
  for( i = 0; i <= n; i++ ) vbmat->bsz[i] = bsz[i];	
  vbmat->wk = (FLOAT *) Malloc(2*nn*sizeof(FLOAT), "setupVBILUT:5" );
  vbmat->lu = NULL;
  vbmat->plu = NULL;
  vbmat->dd1 = NULL;
  vbmat->dd2 = NULL;
  //vbmat->L = (csptr) Malloc(sizeof(SparMat), "setupILUT:6" );
  //if (setupCS(vbmat->L, len,1)) return 1;
  //vbmat->U = (csptr) Malloc(sizeof(SparMat), "setupILUT:7" );
  //if (setupCS(vbmat->U, len,1)) return 1;
  return 0;
}


void vbmatvec1(vbsptr vbmat, FLOAT *x, FLOAT *y )
{
/*-------------------- matrix -- vector product in VB format, for rectangle only
------------------------------y = A*x---------------------------------------*/
  int i, j, nzcount, col, inc = 1, dim, sz, nBs, nBsj; 
  int n = vbmat->n, *ja, *bsz = vbmat->bsz,*bszc = vbmat->bszc;
  FLOAT one=1.0;
  BData *ba;
  
  for( i = 0; i < n; i++ ) {
    nBs = bsz[i];
    dim = B_DIM(bsz,i);
    for( j = 0; j < dim; j++ ) 
      y[nBs+j] = 0;
    nzcount = vbmat->nzcount[i];
    ja = vbmat->ja[i];
    ba = vbmat->ba[i];
    for( j = 0; j < nzcount; j++ ) {
      col = ja[j];
      nBsj = bszc[col];
      sz = B_DIM(bszc,col);
/*-------------------- operation:  y = Block*x + y */
      GGEMV ("n", dim, sz,one, ba[j],dim,&x[nBsj],inc,one,&y[nBs],inc);
    }
  }
}

void iilutptrsolver(iilutptr ilusch, FLOAT *y) 
{
/*---------------------------------------------------------------------
|  Forward solve for Schur complement part = 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return
| y       = solution of LU x = y. [overwritten] 
|---------------------------------------------------------------------*/
//typedef struct IILUTfac {//for coarse point ilutp struct
    //int n;        
    //int *bsz;                  
    //ilutptr *llu; 
//} IIluSpar, *iilutptr;
/*-------------------- local variables                        */
  int n = ilusch->n, j, first=0;//, *perm = ilusch->rperm;
  ilutptr lup = NULL;
  //FLOAT *work = ilusch->wk; 
/*-------------------- begin: right scaling                          */
   //if (ilusch->D1 != NULL) 
     //dscale(n, ilusch->D1, y, y); 
/*-------------------- ONE SIDED ROW PERMS */
   //if (perm != NULL) { 
     for (j=0; j<n; j++){
       //work[perm[j]] = y[j];
	lup = ilusch->llu[j]; 
	//first = bsz[j];
	SchLsol(lup, &y[first]);//&x[first] 
        SchUsol(lup, &y[first]);
	first = first + lup->n;
/*--------------------  L solve proper */
     //Lsol(ilusch->L, work, y); 
     }
   //} else 
     //Lsol(ilusch->L, y, y); 
/*---------------end of iilutptrsolver---------------------------------------
----------------------------------------------------------------------*/
}

void vbmatvecz(vbsptr vbmat, FLOAT *x, FLOAT *y, FLOAT *z )
{
/*---------------------------------------------------------------------
| This function does the matrix vector  z = y - A x.
|----------------------------------------------------------------------
| on entry:
| vbmat  = the matrix, has to be rectangle (in VBSpaFmt form)
| x, y   = two input vector
|
| on return
| z    = the result:  y - A * x
| z-location must be different from that of x 
| i.e., y and x are used but not modified.
|--------------------------------------------------------------------*/
  int i, j, nzcount, col, inc = 1, dimR, dimC, nBs, nBsj; 
  int n = vbmat->n, *ja, *bsz = vbmat->bsz, *bszc = vbmat->bszc;
  FLOAT one = 1.0, minusone = -1.0, *tt;
  BData *ba;
	//printf("\n MAX_BLOCK_SIZE=%d",MAX_BLOCK_SIZE); getchar();
  tt = (FLOAT *) Malloc((MAX_BLOCK_SIZE)*sizeof(FLOAT), "vbmatvecz" );	
  
  for( i = 0; i < n; i++ ) {
    nBs = bsz[i];
    dimR = B_DIM(bsz,i);
    //for( j = 0; j < dim; j++ ) 
      //y[nBs+j] = 0;
    nzcount = vbmat->nzcount[i];
    ja = vbmat->ja[i];
    ba = vbmat->ba[i];
    copyBData( 1, dimR, tt, &y[nBs], 0 );//t = y[i], 

    for( j = 0; j < nzcount; j++ ) {
      col = ja[j];
      nBsj = bszc[col];
      dimC = B_DIM(bszc,col);
/*-------------------- operation:  y = Block*x + y */
      GGEMV ("n", dimR, dimC, minusone, ba[j],dimR,&x[nBsj],inc,one,tt,inc);//t -= kr[k] * x[ki[k]];
    }
    copyBData( 1, dimR, &z[nBs],tt , 0 );
  }
	free(tt);
}


int vbascend (vbp4ptr levmat, FLOAT *x, FLOAT *wk) 
{
/*---------------------------------------------------------------------
| This function does the (block) backward substitution: 
|
|     |            |  |     |    |    |
|     | U  L^{-1}F |  | wk1 |    | x1 |
|     |            |  |     | =  |    |
|     | 0       S  |  | wk2 |    | x2 |  <<-- x2 already computed.
|     |            |  |     |    |    |       and we need x1
|
|    with x2 = S^{-1} wk2 [assumed to have been computed ] 
|--------------------------------------------------------------------*/
 /*--------------------  local variables  */
  int j, len=levmat->bsz[levmat->n], lenB=levmat->bsz[levmat->nB], *qperm=levmat->perm;
  FLOAT *work = levmat->wk; 
  /*-------------------- copy x onto wk */  
  vbmatvec1(levmat->F, &x[lenB], work);   /*  work = F * x_2   */

     if (levmat->lu) 
       vblusolC(work, work,levmat->lu);//ierr = vblusolC(work,wk,lu);
     else
       iilutptrsolver(levmat->plu, work);

  //ierr = vblusolC(work, work,levmat->lu);         /*  work = L \ work   L x = b; b,x */
  for (j=0; j<lenB; j++)               /*  wk1 = wk1 - work  */
    work[j] = x[j] - work[j];
  //Usol(levmat->U, work, work);         /*  wk1 = U \ wk1 */ 
  memcpy(&work[lenB],&x[lenB],(len-lenB)*sizeof(FLOAT));
/*---------------------------------------
|   apply reverse permutation
|--------------------------------------*/
  for (j=0; j<len; j++)
     wk[j] = work[qperm[j]];     
  return 0;
}
/*----end-of-vbascend----------------------------------------------------
|----------------------------------------------------------------------
|--------------------------------------------------------------------*/


int vbdescend(vbp4ptr levmat, FLOAT *x, FLOAT *wk)
{
/*---------------------------------------------------------------------
| This function does the (block) forward elimination in ARMS
|                       new       old
|     |            |  |     |    |    |
|     | L        0 |  | wx1 |    | x1 |
|     |            |  |     | =  |    | 
|     | EU^{-1}  I |  | wx2 |    | x2 |
|     |            |  |     |    |    |
| x used and not touched -- or can be the same as wk.
|--------------------------------------------------------------------*/
/*  local variables   */
  int j, len=levmat->bsz[levmat->n], lenB=levmat->bsz[levmat->nB], *iperm=levmat->rperm; 
  FLOAT *work = levmat->wk; 
  vbiluptr lu = levmat->lu;
  iilutptr plu = levmat->plu;
/*------------------------------------------------------
|   apply permutation P to rhs wk
|-----------------------------------------------------*/
  for (j=0; j<len; j++)
    work[iperm[j]] = x[j] ;
    //output_dblvector("worksecondpart.coo" ,x,0, len);
    //output_dblvector("workpermuted.coo" ,work,0, len);
     //printf("\n levmat->lu=%d,levmat->plu=%d",levmat->lu,levmat->plu); getchar();
//output_dblvector("workbeforep.coo" ,work,0, lenB);	
//output_intvector("bszinB.coo" ,plu->bsz,0, plu->n+1);
     if (lu) 
       vblusolC(work,wk,lu);
     else
       {memcpy(wk,work,len*sizeof(FLOAT));iilutptrsolver(plu,wk);}//void iilutptrsolver(iilutptr ilusch, FLOAT *y) 
//output_dblvector("workafterp.coo" ,wk,0, lenB);
    //ierr = vblusolC(work,wk,lu);//int vblusolC( FLOAT *y, FLOAT *x, vbiluptr lu)
  //Lsol(levmat->L, work, wk);       /* sol:   L x = x                 *///void Lsol(csptr mata, FLOAT *b, FLOAT *x)  L x = b
  //Usol(levmat->U, wk, work);        /* sol:   U work(2) = work         */
/*-------------------- compute x[lenb:.] = x [lenb:.] - E * work(1) */
  //printf("\n lenB=%d",lenB); getchar();
	//output_dblvector("work1.coo" ,work,0, lenB);
	//output_dblvector("work2.coo" ,&work[lenB],0, len-lenB);
	//outputvbmat1(levmat->E,"levmat->E.coo",1);
  vbmatvecz (levmat->E, wk, &work[lenB], &wk[lenB]); //void vbmatvecz(vbsptr vbmat, FLOAT *x, FLOAT *y, FLOAT *z )//z = y - A x
	//output_dblvector("wkde.coo" ,wk,0, len);
  return 0;
}
/*----end-of-vbdescend---------------------------------------------------
|----------------------------------------------------------------------
|--------------------------------------------------------------------*/

//int permuten(int *iwork, int *iworkn, int *bsz, int n)
//{
/*----------------------------------------------------------------------
| Generate the full length permutation array.
|----------------------------------------------------------------------
| on entry:
|==========
|     iwork = block permutation array.
|     n   =  size of iwork
|     bsz = block index of current level
|
| On return:
|===========
|
|
|  iworkn = full length permutation array. 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
 /* int k,i,j;


  for( k=0,i = 0; i < n; i++ ){ 	
  	for(j = 0; j < (bsz[iwork[i]+1]-bsz[iwork[i]]); j++ ){
		iworkn[k] = bsz[iwork[i]]+j;k++;
	}
  }
  return 0;
}*/

void vbdscale(int nn, double *dd, FLOAT *x, FLOAT * y)
{ 
/* Computes  y == DD * x                               */
/* scales the vector x by the diagonal dd - output in y */


  int k;//dimsz,first,i;
  for (k=0; k<nn; k++){
	y[k] = dd[k]*x[k];
	//DDOT(n,x,incx,y,incy)//y[k] = dd[k]*x[k];//dscale(nloc, levmat->D2, &x[first], &x[first]) ;
  }
}

void VBSchLUsol(vbilutptr ilusch, FLOAT *y) 
{
/*---------------------------------------------------------------------
|  Forward solve for Schur complement part = 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return
| y       = solution of LU x = y. [overwritten] 
|---------------------------------------------------------------------*/
/*-------------------- local variables                        */
  int nn = ilusch->bsz[ilusch->n];//j, *perm = ilusch->rperm,
  //FLOAT *work = ilusch->wk; 
/*-------------------- begin: right scaling                          */
//printf("\n ilusch->D1=%d",ilusch->D1); getchar();
   if (ilusch->dd1 != NULL) 
     vbdscale(nn, ilusch->dd1, y, y);

   if (ilusch->D1 != NULL) 
     vbdscale(nn, ilusch->D1, y, y); 
/*-------------------- ONE SIDED ROW PERMS */
   //if (perm != NULL) { 
     //for (j=0; j<n; j++)
       //work[perm[j]] = y[j]; 
/*--------------------  L solve proper */
     //Lsol(ilusch->L, work, y); 
   //} else 
     vblusolC(y, y,ilusch->lu);//Lsol(ilusch->lu, y, y); 
/*-------------------- case when diagonal scaling is done on columns    */ 
    if (ilusch->D2 !=NULL) 
      vbdscale(nn, ilusch->D2, y, y);

    if (ilusch->dd2 !=NULL) 
      vbdscale(nn, ilusch->dd2, y, y);
/*---------------end of SchLsol---------------------------------------
----------------------------------------------------------------------*/
}

void pSchLUsol(vbilutptr ilusch, FLOAT *y) 
{
/*---------------------------------------------------------------------
|  Forward solve for Schur complement part = 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return
| y       = solution of LU x = y. [overwritten] 
|---------------------------------------------------------------------*/
/*-------------------- local variables                        */
  int  nn = ilusch->bsz[ilusch->n];//j, *perm = ilusch->rperm,ierr = 0,
  //FLOAT *work = ilusch->wk; 
/*-------------------- begin: right scaling                          */
//printf("\n ilusch->D1=%d",ilusch->D1); getchar();
   if (ilusch->D1 != NULL) 
     vbdscale(nn, ilusch->D1, y, y); 
/*-------------------- ONE SIDED ROW PERMS */
   //if (perm != NULL) { 
     //for (j=0; j<n; j++)
       //work[perm[j]] = y[j]; 
/*--------------------  L solve proper */
     //Lsol(ilusch->L, work, y); 
   //} else 
    SchLsol(ilusch->plu,y);SchUsol(ilusch->plu,y);
     //ierr = vblusolC(y, y,ilusch->lu);//Lsol(ilusch->lu, y, y); 
/*-------------------- case when diagonal scaling is done on columns    */ 
    if (ilusch->D2 !=NULL) 
      vbdscale(nn, ilusch->D2, y, y);
/*---------------end of SchLsol---------------------------------------
----------------------------------------------------------------------*/
}

int cleanVBILUT(vbilutptr amat, int indic)
{
/*----------------------------------------------------------------------
| Free up memory allocated for IluSpar structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a IluSpar struct.
|  indic    = indicator for number of levels.  indic=0 -> zero level.
|--------------------------------------------------------------------*/

  /*----------------*/
   
  if (amat->wk) {
    free(amat->wk); 
    amat->wk = NULL;
  }
  if (amat->lu) {
    cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    amat->lu = NULL;
  }  //cleanCS(amat->L);
  //cleanCS(amat->U);
  if (amat->plu) {
    cleanILUT(amat->plu,0);//extern int cleanILUT(ilutptr amat, int indic);//cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    free(amat->plu);
    amat->plu = NULL;
  }  //cleanCS(amat->L);
  //cleanCS(amat->U);

  if (indic) cleanVBMat(amat->C);//cleanCS(amat->C);  
/*-------------------- nonsymmetric permutation */
  if (amat->rperm) {
    free(amat->rperm);
    amat->rperm = NULL;
  }
  if (amat->perm) {
    free(amat->perm); 
    amat->perm = NULL;
  }
  
/*-------------------- ilutp permutation */
 // if (amat->perm2) free(amat->perm2);
/*-------------------- diagonal scalings */
    /*if( amat->D1 ) {
        for( i = 0; i < amat->n; i++ ) {
            if( amat->D1[i] ) free( amat->D1[i] );
        }
        free( amat->D1 );
    }
    if( amat->D2 ) {
        for( i = 0; i < amat->n; i++ ) {
            if( amat->D2[i] ) free( amat->D2[i] );
        }
        free( amat->D2 );
    }*/
//~ printf("\n amat->D1=%d,amat->D2=%d",amat->D1,amat->D2);
    if (amat->D1) {free(amat->D1);amat->D1 = NULL;}
    if (amat->D2) {free(amat->D2);amat->D2 = NULL;}
    if (amat->dd1) {free(amat->dd1);amat->dd1 = NULL;}
    if (amat->dd2) {free(amat->dd2);amat->dd2 = NULL;}
    if( amat->bsz ) {free( amat->bsz );amat->bsz = NULL;}
   return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBILUT
|--------------------------------------------------------------------*/

int cleanVBP4(vbp4ptr amat)
{
/*----------------------------------------------------------------------
| Free up memory allocated for VBPer4Mat structs.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a VBPer4Mat struct.
|--------------------------------------------------------------------*/

/*  -------------------------- */


  if (amat == NULL) return 0;
  if (amat->n < 1) return 0;
   

  if (amat->perm) {
    if (amat->perm) free(amat->perm); 
    amat->perm = NULL;
  }
  
  if (!amat->symperm) { 
    if (amat->rperm) free(amat->rperm); 
    amat->rperm = NULL;
  } 
  
  if (amat->F) {
    cleanVBMat1(amat->F); //cleanVBMat(vbsptr vbmat);
    amat->F = NULL;
  }
  if (amat->E) {
    cleanVBMat1(amat->E); 
    amat->E = NULL;
  }
  if (amat->lu) {
    cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    amat->lu = NULL;
  }

  if (amat->plu) {
    cleanIILUT(amat->plu);//extern int cleanILUT(ilutptr amat, int indic);//cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    amat->plu = NULL;
  } 
  //if (amat->plu) {
    //cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    //amat->lu = NULL;
  //}
  //if (amat->L) {
    //cleanCS(amat->L);
    //amat->L = NULL;
   //}
  //if (amat->U) {
    //cleanCS(amat->U);
    //amat->U = NULL;
  //}
  
  if (amat->prev == NULL) 
    if (amat->wk) free(amat->wk);  
  
    /*if( amat->D1 ) {
        for( i = 0; i < amat->n; i++ ) {
            if( amat->D1[i] ) free( amat->D1[i] );
        }
        free( amat->D1 );
    }
    if( amat->D2 ) {
        for( i = 0; i < amat->n; i++ ) {
            if( amat->D2[i] ) free( amat->D2[i] );
        }
        free( amat->D2 );
    }*/
    if (amat->D1) {free(amat->D1);amat->D1 = NULL;}
    if (amat->D2) {free(amat->D2);amat->D2 = NULL;}
    if( amat->bsz ) {free( amat->bsz );amat->bsz = NULL;}
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanVBP4
|--------------------------------------------------------------------*/


int cleanBARMS(vbarms ArmsPre)
{
  vbp4ptr amat = ArmsPre->levmat;
  vbilutptr cmat = ArmsPre->ilus;
/*----------------------------------------------------------------------
| Free up memory allocated for entire BARMS preconditioner.
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )  =  Pointer to a Per4Mat struct.
| ( cmat )  =  Pointer to a IluSpar struct.
|--------------------------------------------------------------------*/
/* case when nlev == 0 */  
  int indic=(amat->nB != 0) ;
    /*  && amat->next !=NULL) ; */
  
  vbp4ptr levc, levn;

  levc = amat; 

  if (indic) { 
    while (levc) {
      if (cleanVBP4(levc)) return(1) ; 
      levn = levc->next;
      free(levc);
      levc = levn;
    }		
  }	
   else 	
     if (amat) {
       free(amat) ; 
       amat = NULL;
     }
  
  cleanVBILUT(cmat,indic); 
  
  
  if (cmat) {
    free(cmat);	
    cmat = NULL;
  }
 // if (ArmsPre) {
    //free(ArmsPre);	
    //ArmsPre = NULL;
  //}
  
  return 0;
}
/*---------------------------------------------------------------------
|     end of cleanBARMS 
|--------------------------------------------------------------------*/

int nnz_iilutptr(iilutptr lu )
{  
  /*-------------------- counts number of nonzero entries in preconditioner */
    //int *bsz = lu->bsz;
    int nnz = 0,nzblock = 0,i;//, i, j, col;nzcount, 
    ilutptr lup;
    for( i = 0; i < lu->n; i++ ) {
	lup = lu->llu[i];

	nzblock = cs_nnz(lup->L)+cs_nnz(lup->U);
        //nzcount = 0;
        //for( j = 0; j < lu->L->nzcount[i]; j++ ) {
            //col = lu->L->ja[i][j];
            //nzcount += B_DIM(bsz,col);
        //}
        //for( j = 0; j < lu->U->nzcount[i]; j++ ) {
            //col = lu->U->ja[i][j];
            //nzcount += B_DIM(bsz,col);
        //}
        //nzcount += B_DIM(bsz,i);  /* diagonal */
        //nzcount *= B_DIM(bsz,i);
        nnz += nzblock;
    }
    return nnz;
}


int nnz_vbilu(vbiluptr lu )
{  
  /*-------------------- counts number of nonzero entries in preconditioner */
    int *bsz = lu->bsz;
    int nzcount, nnz = 0, i, j, col;
    for( i = 0; i < lu->n; i++ ) {
        nzcount = 0;
        for( j = 0; j < lu->L->nzcount[i]; j++ ) {
            col = lu->L->ja[i][j];
            nzcount += B_DIM(bsz,col);
        }
        for( j = 0; j < lu->U->nzcount[i]; j++ ) {
            col = lu->U->ja[i][j];
            nzcount += B_DIM(bsz,col);
        }
        nzcount += B_DIM(bsz,i);  /* diagonal */
        nzcount *= B_DIM(bsz,i);
        nnz += nzcount;
    }
    return nnz;
}

int vblev4_nnz(vbp4ptr levmat, int *lev, FILE *ft) 
{
  /* counts all nonzero elements in levmat struct  -- 
     recursive */
  int nnzT, nnzLDU, nnzF, nnzE, nnzDown=0;//, nnzU
  vbp4ptr nextmat;

//~ printf("levmat->lu = %p",levmat->lu);
     //if (levmat->lu) 
       nnzLDU = nnz_vbilu(levmat->lu );
     //else
       //nnzLDU = nnz_iilutptr(levmat->plu );//int nnz_iilutptr(iilutptr lu );//{nnzLDU = cs_nnz(levmat->plu->L)+cs_nnz(levmat->plu->U);} 

  //nnzLDU = nnz_vbilu(levmat->lu );//nnzL = cs_nnz(levmat->L); int nnz_vbilu(vbiluptr lu )
  //nnzU = cs_nnz(levmat->U); 
  nnzF = memVBMat1(levmat->F);//int memVBMat1( vbsptr vbmat ) 
  nnzE = memVBMat1(levmat->E); 
  nnzT = nnzLDU+nnzF+nnzE;
  /* print */ 
#if 0
  if (*lev == 0)
    fprintf(ft,
  	    "\nnnz/lev used:      LU        F        E    subtot\n");
  fprintf(ft,"    Level %2d %8d %8d %8d %8d\n",
  	  *lev, nnzLDU, nnzF, nnzE, nnzT);
#endif
  if (*lev == 0)
    printf("\nnnz/lev used:      LU        F        E    subtot\n");
  printf("    Level %2d %8d %8d %8d %8d\n",
  	  *lev, nnzLDU, nnzF, nnzE, nnzT);
  (*lev)++;
  nextmat = levmat->next; 
  if (nextmat != NULL) 
   nnzDown = vblev4_nnz(nextmat, lev, ft);
   return (nnzT+nnzDown); 
}

int nnz_vbarms (vbarms PreSt,  FILE *ft)
{ 
/*-------------------------------------------------------
| computes and prints out total number of nonzero elements
| used in BARMS factorization 
+--------------------------------------------------------*/
  vbp4ptr levmat   = PreSt->levmat; 
  vbilutptr ilschu = PreSt->ilus; 
  int nlev       = PreSt->nlev;
  int ilev=0,nnz_lev,nnz_sch,nnz_tot; 
  nnz_lev = 0; 
  if (nlev) nnz_lev+= vblev4_nnz(levmat, &ilev, ft);

     //~ if (levmat->lu) 
       nnz_sch = nnz_vbilu(ilschu->lu );
     //~ else
       //~ {nnz_sch = cs_nnz(ilschu->plu->L)+cs_nnz(ilschu->plu->U);} 

  //nnz_sch = nnz_vbilu(ilschu->lu );//nnz_sch = cs_nnz(ilschu->L)+cs_nnz(ilschu->U);
//printf("\n nnz_sch=%d",nnz_sch); //getchar();	//outputvbmat(ilschu->L,"lastBL.coo",1);//outputvbmat(schur,"schur.coo",1);
	//outputvbmat(ilschu->U,"lastBU.coo",1);
	//outputvbD( lu, "DD.coo", 1);
  if (nlev) nnz_sch += memVBMat(ilschu->C);
  nnz_tot = nnz_lev+nnz_sch; 
#if 0
  fprintf(ft,"\n");
  fprintf(ft,"Total nonzeros for interm. blocks.... =  %10d\n",nnz_lev);
  fprintf(ft,"Total nonzeros for last level ....... =  %10d\n",nnz_sch);
  fprintf(ft,"Grand total.......................... =  %10d\n",nnz_tot);
#endif 
  printf("\n");
  printf("Total nonzeros for interm. blocks.... =  %10d\n",nnz_lev);
  printf("Total nonzeros for last level ....... =  %10d\n",nnz_sch);
  printf("Grand total.......................... =  %10d\n",nnz_tot);
  return nnz_tot;
}

void vproduct(int n, FLOAT *dst, FLOAT *src, int isig)
{ 
  int i;
  if( isig == 0 )
    for( i = 0; i < n; i++ ) dst[i] = src[i]*dst[i];
  else
    for( i = 0; i < n; i++ ) dst[i] = -src[i]*dst[i];
}

int vbfroscalC(vbsptr mata, double *diag, int nrm)
{
    /*---------------------------------------------------------------------
|
| This routine scales each row of mata so that the norm is 1, do the point operation on the blocked matrix.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in VBSpaFmt form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> row j is a zero row
|--------------------------------------------------------------------*/
    /*   local variables    */
    int i, j, k, m, col=0 ,dimR, dimC, blocksz, first, first1, rindex;
    int *bsz = mata->bsz,*ja;
    int max_row_sz = MAX_BLOCK_SIZE*sizeof(double);
    double *scal;//t;
    BData *ba;

    scal = (double*)Malloc( max_row_sz, "vbroscalC" );

    for (i=0; i<mata->n; i++) {
        memset(scal, 0, max_row_sz);//scal = 0.0;//memset(mat[i], 0, nnr*sizeof(FLOAT));
        ba = mata->ba[i];
        ja = mata->ja[i];
        dimR = B_DIM(bsz,i);
        if (nrm == 0) {
            for (j=0; j<mata->nzcount[i]; j++){
                col = ja[j];
                dimC = B_DIM(bsz,col);
                blocksz = dimR*dimC;
                //t = vbnorm2(blocksz,ba[j]);
                //scal = max(t,scal);//if (t > scal) scal = t;//fabs(ba[j]) > scal
                for (k=0; k<blocksz; k++){
                    rindex = k%dimR;
                    if (ABS_VALUE(ba[j][k]) > scal[rindex])
                        scal[rindex] = ABS_VALUE(ba[j][k]);//if (ABS_VALUE(kr[k]) > scal) scal = ABS_VALUE(kr[k]);k%dimR//k/dimR
                }
            }
        }
        else if (nrm == 1) {
            for (j=0; j<mata->nzcount[i]; j++){
                col = ja[j];
                dimC = B_DIM(bsz,col);
                blocksz = dimR*dimC;
                for (k=0; k<blocksz; k++){
                    rindex = k%dimR;
                    scal[rindex] += ABS_VALUE(ba[j][k]);
                }
                //t = vbnorm2( blocksz,ba[j]);
                //scal += t;//scal += fabs(ba[j]);
            }
        }
        else {  /* nrm = 2 */
            for (j=0; j<mata->nzcount[i]; j++){
                col = ja[j];
                dimC = B_DIM(bsz,col);
                blocksz = dimR*dimC;
                for (k=0; k<blocksz;k++){
                    rindex = k%dimR;
                    scal[rindex] += ba[j][k]*ba[j][k];
                }
                //t = vbnorm2( blocksz,ba[j]);
                //scal += t*t;//scal += ba[j]*ba[j];
            }
        }
        if (nrm == 2) for (j=0; j<dimR;j++) scal[j] = sqrt(scal[j]);
        for (j=0; j<dimR;j++) {
            if (scal[j] == 0.0)
                scal[j] = 1.0;
                /* YS. return i+1; */

            else
                scal[j] = 1.0 / scal[j];
        }
        first = bsz[i];
        memcpy(&diag[first], scal, dimR*sizeof(double));
        //for (j=0; j<dimR; j++){
        //diag[first+j] = scal;//memcpy(&diag[first],w,dimR*sizeof(FLOAT));
        //}
        /*void generatefull(FLOAT *d1, FLOAT *d11, int n, int *bsz)
{//subroutine for testing 
int i,j,dimsz,first;
    for( i = 0; i < n; i++ ) {
        dimsz = B_DIM(bsz,i);
        first = bsz[i];
        for(j=0; j<dimsz;j++)
        d11[first+j] = d1[i];
        //printf("\n rhs[i]=%20.16e,rhs1[iworkn[i]]=%20.16e,iworkn[i]=%d,i=%d",rhs[i],rhs1[iworkn[i]],iworkn[i],i); getchar();
    }
}*/
        //printf("\n i=%d,scal=%f",i,scal); getchar();
        for (j = 0; j < mata->nzcount[i]; j++){
            col = ja[j];
            dimC = B_DIM(bsz,col);
            blocksz = dimR*dimC;
            for (k=0; k<dimC;k++){
                first1 = k*dimR;
//                output_dblvector("&ba[j][first1]before.coo" ,&ba[j][first1] ,0 , dimR);
                //output_dblvector("&scalbefore.coo" ,scal ,0 , dimR);
                for (m = 0; m < dimR; m++)
                     ba[j][first1+m] = ba[j][first1+m]*scal[m];
//                vproduct(dimR, &ba[j][first1], scal, inc);//void vproduct(int n, FLOAT *dst, FLOAT *src, int isig)//DDOT(dimR,&ba[j][first1],inc,scal,inc);//DDOT(n,x,incx,y,incy)
//                output_dblvector("&ba[j][first1]after.coo" ,&ba[j][first1] ,0 , dimR);
                //output_dblvector("&scalafter.coo" ,scal ,0 , dimR);
                //rindex = k%dimR;
                //scal[rindex] += ba[j][k]*ba[j][k];
            }
            //printf("\n i=%d,scal=%f,dimC = %d,blocksz=%d",i,scal,dimC,blocksz); getchar();
            //output_dblvector("ba[j]before.coo" ,ba[j],0, blocksz);
            //DSCAL(blocksz,scal,ba[j],inc);//DSCAL(blocksz,scal,ba[k],inc);//ba[k] = ba[k] * scal;//DSCAL(n,alpha,x,incx)???
            //printf("\n i=%d,scal=%f,dimC = %d,blocksz=%d",i,scal,dimC,blocksz); getchar();
            //output_dblvector("ba[k]after.coo" ,ba[k],0, blocksz);
            /*for (k=0; k<dimR;k++){
         //first1 = k*dimR;
        inc = blocksz-k;
         DSCAL(inc,scal[k],&ba[j][k],dimC);//vproduct(dimR,&ba[j][first1],scal,inc);//void vproduct(int n, FLOAT *dst, FLOAT *src, int isig)//DDOT(dimR,&ba[j][first1],inc,scal,inc);
         //rindex = k%dimR;
         //scal[rindex] += ba[j][k]*ba[j][k];
     }*/
        }
    }
    free(scal);
    return 0;
}
/*---------------end of vbfroscalC-----------------------------------------
----------------------------------------------------------------------*/


int vbfcoscalC(vbsptr mata, double *diag, int nrm)
{
    /*---------------------------------------------------------------------
|
| This routine scales each column of mata so that the norm is 1, do the point operation on the blocked matrix.
|
|----------------------------------------------------------------------
| on entry:
| mata  = the matrix (in VBSpaFmt form)
| nrm   = type of norm
|          0 (\infty),  1 or 2
|
| on return
| diag  = diag[j] = 1/norm(row[j])
|
|     0 --> normal return
|     j --> column j is a zero column
|--------------------------------------------------------------------*/
    /*   local variables    */
    int i, col, k, j, dimR, dimC, inc = 1,blocksz,nn = mata->bsz[mata->n],first1;
    BData *kr;
    int *bsz = mata->bsz, *ki;

#if defined(DBL_CMPLX)
    FLOAT diagEntry = 0;
#endif
    for (i=0; i<nn; i++)
        diag[i] = 0.0;
    /*---------------------------------------
|   compute the norm of each column
|--------------------------------------*/
    for (i=0; i<mata->n; i++) {
        kr = mata->ba[i];
        ki = mata->ja[i];
        dimR = B_DIM(bsz,i);
        if (nrm == 0) {
            for (k=0; k<mata->nzcount[i]; k++) {
                /*col = ja[k];
        dimC = B_DIM(bsz,col);
        blocksz = dimR*dimC;*/
                col = ki[k];
                dimC = B_DIM(bsz,col);
                blocksz = dimR*dimC;
                for (j=0; j<blocksz; j++)
                    diag[bsz[col]+j/dimR] = max(ABS_VALUE(kr[k][j]), diag[bsz[col]+j/dimR]);
                //t = vbnorm2(blocksz,kr[k]);//bsz[col]+k/dimR
                //diag[col] = max(t,diag[col]);//if (t > diag[col]) diag[col] = fabs(kr[k]);//if (fabs(kr[k]) > diag[col]) diag[col] = fabs(kr[k]);
            }
        }
        else if (nrm == 1) {
            for (k=0; k<mata->nzcount[i]; k++){
                col = ki[k];
                dimC = B_DIM(bsz,col);
                blocksz = dimR*dimC;
                for (j=0; j<blocksz; j++)
                    diag[bsz[col]+j/dimR] += ABS_VALUE(kr[k][j]);
                //t = vbnorm2(blocksz,kr[k]);
                //diag[col] += t;//diag[ki[k]] += fabs(kr[k]);
            }
        }
        else {  /*  nrm = 2 */
            for (k=0; k<mata->nzcount[i]; k++){
                col = ki[k];
                dimC = B_DIM(bsz,col);
                blocksz = dimR*dimC;
                for (j=0; j<blocksz; j++)
                    diag[bsz[col]+j/dimR] += kr[k][j]*kr[k][j];
                //t = vbnorm2(blocksz,kr[k]);
                //diag[col] += t*t;//diag[ki[k]] += kr[k]*kr[k];
            }
        }
    }
    if (nrm == 2) {
        for (i=0; i<nn; i++)
            diag[i] = sqrt(diag[i]);
    }
    /*---------------------------------------
|   invert
|--------------------------------------*/
    for (i=0; i<nn; i++) {
        if (diag[i] == 0.0)
            /* return i+1;*/
            diag[i] = 1.0;
        else
            diag[i] = 1.0 / diag[i];
    }
    /*---------------------------------------
|   C = A * D
|--------------------------------------*/
    for (i=0; i<mata->n; i++) {
        kr = mata->ba[i];
        ki = mata->ja[i];
        dimR = B_DIM(bsz,i);
        //first = bsz[i];
        /*dimsz = B_DIM(bsz,i);
        first = bsz[i];
        for(j=0; j<dimsz;j++)
        d11[first+j] = d1[i];
        //printf("\n rhs[i]=%20.16e,rhs1[iworkn[i]]=%20.16e,iworkn[i]=%d,i=%d",rhs[i],rhs1[iworkn[i]],iworkn[i],i); getchar();
    }*/
        for (k=0; k<mata->nzcount[i]; k++){
            col = ki[k];
            dimC = B_DIM(bsz,col);
            //blocksz = dimR*dimC;
            //printf("\n i=%d,scal=%f,dimC = %d,blocksz=%d",i,scal,dimC,blocksz); getchar();
            for (j=0; j<dimC; j++){
                first1 = j*dimR;
#if defined(DBL_CMPLX)
                diagEntry = (FLOAT)diag[bsz[col]+j];
                GSCAL(dimR, diagEntry, &kr[k][first1], inc);//kr[k] = kr[k] * diag[ki[k]];
#else
                GSCAL(dimR, diag[bsz[col]+j],&kr[k][first1],inc);//kr[k] = kr[k] * diag[ki[k]];
#endif

            }
        }
    }
    return 0;
}


void luinv( int n, FLOAT *a, FLOAT *x, FLOAT *y )
{
/*--------------------------------------------------------
 *    does the operation y = inv(a) * x
 *    where a has already been factored by Gauss.
 *    LUy = x
 *------------------------------------------------------*/
    int i, j, bsA, bsB;
    FLOAT sum;
    /* Ly0 = x -- use Lsol ? */   
    for( i = 0; i < n; i++ ) {
        sum = x[i];
        bsA = i - n;
        for( j = 0; j < i; j++ ) {
            bsA += n;
            sum -= a[bsA] * y[j]; /* a(i,j) * y(j) */
        }
        y[i] = sum;
    }
    /* Uy = y0 */
    bsB = i * n;
    for( i = n-1; i >= 0; i-- ) {
        sum = y[i];
        bsB -= n;
        bsA = i+bsB;
        for( j = i+1; j < n; j++ ) {
            bsA += n;
            sum -= a[bsA] * y[j]; /* a(i,j) * y(j) */
        }
        y[i] = sum * a[bsB+i]; /* a(i,i) */
    }
}

int indsetC2(csptr mat, int bsize, int *iord, int *nnod, double tol,
	    int nbnd, int *nBB, int *nset) 
{
/*--------------------------------------------------------------------- 
| greedy algorithm for independent set ordering -- 
|----------------------------------------------------------------------
|     Input parameters:
|     -----------------
|     (mat)  =  matrix in SparRow format
|     
|     bsize  =  integer (input) the target size of each block.
|               each block is of size >= bsize. 
|
|     w      =  weight factors for the selection of the elements in the
|               independent set. If w(i) is small i will be left for the
|               vertex cover set. 
|
|     tol    =  a tolerance for excluding a row from independent set.
|
|     nbnd   =  number of interior variables.
|
|     Output parameters:
|     ------------------ 
|     iord   = permutation array corresponding to the independent set 
|     ordering.  Row number i will become row number iord[i] in 
|     permuted matrix.
|     
|     nnod   = (output) number of elements in the independent set.
|     nBB = the array stores the size of independent set.
|     nset = the length of nBB.
|     
|----------------------------------------------------------------------- 
|     the algorithm searches nodes in lexicographic order and groups
|     the (BSIZE-1) nearest nodes of the current to form a block of
|     size BSIZE. The current algorithm does not use values of the matrix.
|---------------------------------------------------------------------*/ 
/*   local variables   */
   int nod, jcount, lastlev, begin, last0, last, nback, mid,
     j1, j2, jcol, inod, jnod, j, k, jcount0, begin0, *rowj, i = 0;
   int prog, n=mat->n, *riord;
   double *w;
   csptr matT,gmat;  

/*-----------------------------------------------------------------------*/
   riord = (int *) Malloc(n*sizeof(int), "indsetC:1" );
   w     = (double *) Malloc(n*sizeof(double), "indsetC:2" );
   matT  = (csptr) Malloc(sizeof(SparMat), "indsetC:3" );
/*  	 call weights to compute the weights for  input matrix.. */
   setupCS(matT, mat->n,1);
   SparTran(mat, matT, 1, 0);
   SparTran(matT, mat, 1, 1); 
   weightsC(mat, w); 
/*---------------------------------------------------------------------- 
| scan all nodes first to eliminate those not satisfying DD criterion 
+----------------------------------------------------------------------*/
   // nbnd = n;//test
   printf("nbnd = %d, in indsetC2 \n",nbnd);

   nback = n-1; 
   nod = 0;
   for(j=0; j<n; j++)
     iord[j] = -1; 
   for(j = n-1; j >=nbnd; j--) {
     add2com(&nback, j, iord, riord);
   }
   for(j=0; j< nbnd; j++) {
     if (w[j] < tol) {
       add2com(&nback, j, iord, riord);
       nod++;
     }
   }
   last = -1;
   for (nod=0; nod<n; nod++) {	
     while (iord[nod] != -1)   {
       if (++nod >= mat->n) goto label50;
     }
/*-------------------- initialize level-set - contains nod (only)*/
     add2is(&last, nod, iord, riord);
     begin   = last;
     begin0  = begin; 
     lastlev = begin;
     jcount  = 1;
/*----------------------------------------------------------------------
|     put all the nearest neighbor nodes of the current node into
|     the block until the number is BSIZE.
|---------------------------------------------------------------------*/
     prog = 1;
     while (jcount < bsize && prog) {
/*--------------------   traverse all the current level-set   */
       last0 = last;
       jcount0 = jcount;
       for (inod=begin; inod<=last0; inod++) {
	 jnod = riord[inod]; 
/*--------------------   This assumes A is not symmetric.   */
	 gmat = mat; 
	 for (k=0; k<2; k++) {
	   rowj = gmat->pj[jnod];
	   for (j=0; j<gmat->nnzrow[jnod]; j++) {
	     jcol = rowj[j];
	     if (iord[jcol] == -1 ) {	
	       add2is(&last, jcol, iord, riord);
	       jcount++;
	     }
	   }
	   gmat = matT; 	
	 }
       }
       prog = jcount > jcount0 ? 1 : 0;
       lastlev = begin;
       begin = last0+1;
     }
     nBB[i] = jcount;
     i++;
/*-----------------------------------------------------------------------
| the neighbors of elements of last level go to the complement   
| gmat loop over original matrix and its transpose 
+-----------------------------------------------------------------------*/ 
     gmat = mat; 
     for (k=0; k<2; k++) {
       for (inod=lastlev; inod<=last; inod++)  {	
	 jnod = riord[inod]; 
	 rowj = gmat->pj[jnod];
	 for (j=0; j<gmat->nnzrow[jnod]; j++){	
	   jcol = rowj[j];
	   if (iord[jcol] == -1) 
	     add2com(&nback, jcol, iord, riord);
	 }
       }
       gmat = matT; 	
     }
/*   reverse ordering for this level   */
     mid = (begin0+last) / 2;
     for (inod=begin0; inod<=mid; inod++) {
       j = last - inod + begin0;
       jnod = riord[inod];
       riord[inod] = riord[j];
       riord[j] = jnod;
     }
   }
/*--------------------------------------------------
|  end-main-loop
|-------------------------------------------------*/
/*-------------------- relabel nodes of vertex cover   */
label50:
   *nnod = last;
   *nset = i;
   j1 = *nnod;
   for (j2=*nnod+1; j2<n; j2++) { 
     if (iord[riord[j2]] > -1) {
       if (++j1 != j2) {
	 j = riord[j2];
	 riord[j2] = riord[j1];
	 riord[j1] = j;
       }
     }
   }
/*-------------------- obtain reverse permutation array   */
   for (j=0; j<n; j++)
     iord[riord[j]] = j;
   (*nnod)++;
   cleanCS(matT); 
   free(riord);
   free(w);
   return 0;
}
/*---------------------------------------------------------------------
|-----end-of-indsetC---------------------------------------------------
|--------------------------------------------------------------------*/



int vbilutD(vbsptr schur, double *droptol, int *lfil, vbilutptr vbmat)
{
/*---------------------------------------------------------------------
| Convert permuted vbspaFmt struct to VBPerMat4 struct 
|                - matrix already permuted
|----------------------------------------------------------------------
| on entry:
|========== 
| ( vbmat )  =  Matrix stored in vbspaFmt format.
|              Internal pointers (and associated memory) destroyed before
|              return.
|
| On return:
|===========
|
| B, E, F, C = 4 blocks in 
| 
|          | B   F |      
|   vbmat= |       | 
|          | E   C | 
| 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
   int i, ierr = 0, max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(FLOAT);
   BData *w;
   vbiluptr lu = NULL;
   //vbsptr FF = NULL,FFF = NULL;
   
/*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/
    w = (BData *)Malloc(schur->n*sizeof(BData),"main");
    for( i = 0; i < schur->n; i++ )
      w[i] = (FLOAT *)Malloc( max_blk_sz, "main" );
    //printf("\n lfil=%d,tol=%20.16e",lfil,tol); getchar();
    lu = (vbiluptr)Malloc( sizeof(VBILUSpar), "main" );
    lu->DiagOpt = 0;
    printf("lfil[7] value is %d\n", lfil[7]);//%f %p %s %c

    if (lfil[7] == -1)
        ierr = vbilutC(schur, lu, lfil[5], droptol[5], w, stderr);//int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp )
    else
        ierr = vbilukC(lfil[7], schur, lu, stderr); //ierr = setupVBP4 (vbmat, schur->n, C->n, F, E , int *bsz,vbiluptr lu );//int setupVBP4 (vbp4ptr vbmat, int Bn, int Cn,  vbsptr F,  vbsptr E , int *bsz,vbiluptr lu )//vbchange

    vbmat->lu = lu;
    
    for( i = 0; i < schur->n; i++ )
      free( w[i] );
    free( w );
    
    if (ierr) return(1); 
   
   // cleanVBMat1( FF );
    //cleanVBMat( FFF ); 	
   return 0;
}


void vbLsolp(int start, vbiluptr lu, FLOAT *y, FLOAT *x)
{
  /*---------------------------------------------------------------------
    |
    | This routine does the forward solve L y = b. where L is upper or
    | bottom part of local low triangular matrix
    |
    | Can be done in place.
    |
    | Zhongze Li, Aug. 17th, 2001
    |
    |----------------------------------------------------------------------
    | on entry:
    | start = the index of the first component
    | n     = one ore than the index owned by the processor 
    | b     = a vector
    | mata  = the matrix (in VBSparMat form)
    |
    | on return
    | y     = the product L^{-1} * b
    |
    |--------------------------------------------------------------------*/
    vbsptr mata = lu->L;
    int n = lu->n, *bsz = lu->bsz, i, j, icol, dim, sz;
    int nzcount, nBs, *ja, inc = 1;
    FLOAT *data, alpha = -1.0, beta = 1.0;
    BData *ba;

    /* Block L solve */
    for( i = 0; i < n; i++ ) {
        dim = B_DIM(bsz, i);
        nBs = bsz[i];
	if (x != y)
	    memcpy(x+nBs, y+nBs, sizeof(FLOAT)*dim);
        nzcount = mata->nzcount[i];
        ja = mata->ja[i];
        ba = mata->ba[i];
        for( j = 0; j < nzcount; j++ )
        {
            icol = ja[j];
            if (bsz[icol] < start) {
                sz = B_DIM(bsz, icol);
                data = ba[j];
                GGEMV("n",  dim,  sz,  alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc);
            }
        }
    }
}

void vbUsolp(int start, vbiluptr lu, FLOAT *y, FLOAT *x)
{
  /*---------------------------------------------------------------------
    |
    | This routine does the backward solve U x = y, where U is upper or
    | bottom part of local upper triangular matrix
    |
    | Can be done in place.
    |
    | Zhongze Li, Aug. 17th, 2001
    |
    |----------------------------------------------------------------------
    | on entry:
    | start = the index of the first component
    | n     = one ore than the index owned by the processor 
    | b     = a vector
    | mata  = the matrix (in VBSparMat form)
    |
    | on return
    | x     = the product U^{-1} * y
    |
    |--------------------------------------------------------------------*/
    vbsptr mata = lu->U;
    int n = lu->n, *bsz = lu->bsz, i, j, icol, dim, sz;
    int nzcount, nBs, *ja, inc = 1, OPT;
    FLOAT *data, alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0;
    BData *ba, *D = lu->D;

    /* Block -- U solve */
    OPT = lu->DiagOpt;
    for( i = n-1; i >= 0; i-- ) {
        nBs = bsz[i];
	if (nBs >= start)
	    continue;

        dim = B_DIM(bsz,i);
	if (x != y)
	    memcpy(x+nBs, y+nBs, sizeof(FLOAT)*dim);
        nzcount = mata->nzcount[i];
        ja = mata->ja[i];
        ba = mata->ba[i];
        for( j = 0; j < nzcount; j++ ) {
            icol = ja[j];
            sz = B_DIM(bsz,icol);
            data = ba[j];
            GGEMV( "n", dim, sz, alpha, data, dim, x+bsz[icol], inc,
		   beta, x+nBs, inc ); 
        }
        data = D[i];
	if (OPT == 1) 
	  luinv( dim, data, x+nBs, lu->bf );
	else
      GGEMV( "n", dim, dim, alpha2, data, dim, x+nBs, inc, beta2, lu->bf, inc);
	
        memcpy(x+nBs, lu->bf, sizeof(FLOAT)*dim);
    }
}

void vbinvsp(int start, vbiluptr lu, FLOAT *y, FLOAT *x)
{
  /*---------------------------------------------------------------------
    |
    | This routine does the backward solve U x = y, where U is upper or
    | bottom part of local upper triangular matrix
    |
    | Can be done in place.
    |
    | Zhongze Li, Aug. 17th, 2001
    |
    |----------------------------------------------------------------------
    | on entry:
    | start = the index of the first component
    | n     = one ore than the index owned by the processor 
    | y     = a vector
    | ilusch  = the LU matrix as provided from the ILU routines.
    |
    | on return
    | x     = the product U^{-1} * y
    |
    |---------------------------------------------------------------------*/
    int n = lu->n, *bsz = lu->bsz, i, j, icol, dim, sz;
    int nzcount, nBs, icolp, *ja, inc = 1, OPT;
    FLOAT *data, alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0;
    vbsptr L, U;
    BData *D, *ba;

    L = lu->L;
    U = lu->U;
    D = lu->D;
    OPT = lu->DiagOpt;
    /* Block L solve */
    for( i = 0; i < n; i++ ) {
        nBs = bsz[i]-start;
	if (nBs < 0)
	    continue;

        dim = B_DIM(bsz,i);
	if (x != y)
	    memcpy(x+nBs, y+nBs, sizeof(FLOAT)*dim);

        nzcount = L->nzcount[i];
        ja = L->ja[i];
        ba = L->ba[i];
        for( j = 0; j < nzcount; j++ ) {
            icol = ja[j];
	    icolp = bsz[icol];
	    if(icolp >= start && icol < n) {
            sz = B_DIM(bsz,icol);
            data = ba[j];
            GGEMV( "n",  dim,  sz,  alpha, data, dim, x+icolp-start,
		   inc, beta, x+nBs, inc );//DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)//y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y
	    }
        }
    }
    /* Block -- U solve */
    for( i = n-1; i >= 0; i-- ) {
        nBs = bsz[i]-start;
	if (nBs < 0)
	    break;

        dim = B_DIM(bsz,i);
        nzcount = U->nzcount[i];
        ja = U->ja[i];
        ba = U->ba[i];
        for( j = 0; j < nzcount; j++ ) {
            icol = ja[j];
            sz = B_DIM(bsz,icol);
            data = ba[j];
            GGEMV( "n", dim, sz, alpha, data, dim, x+bsz[icol]-start, inc,
		   beta, x+nBs, inc ); 
        }
        data = D[i];
	if (OPT == 1) 
	  luinv( dim, data, x+nBs, lu->bf );
	else
      GGEMV( "n", dim, dim, alpha2, data, dim, x+nBs, inc, beta2,
         lu->bf, inc ); //GGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)//y := alpha*A*x + beta*y, or y := alpha*A**T*x + beta*y
	
        memcpy(x+nBs, lu->bf, sizeof(FLOAT)*dim);
    }
}
void VBSchLsol(vbilutptr ilusch, FLOAT *y) 
{
/*---------------------------------------------------------------------
|  Forward solve for Schur complement part = 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return
| y       = solution of LU x = y. [overwritten] 
|---------------------------------------------------------------------*/
/*-------------------- local variables                        */
  int nn = ilusch->bsz[ilusch->n];//j, *perm = ilusch->rperm,
  //FLOAT *work = ilusch->wk; 
/*-------------------- begin: right scaling                          */
//printf("\n ilusch->D1=%d",ilusch->D1); getchar();
   if (ilusch->dd1 != NULL) 
     vbdscale(nn, ilusch->dd1, y, y);

   if (ilusch->D1 != NULL) 
     vbdscale(nn, ilusch->D1, y, y); 
/*-------------------- ONE SIDED ROW PERMS */
   //if (perm != NULL) { 
     //for (j=0; j<n; j++)
       //work[perm[j]] = y[j]; 
/*--------------------  L solve proper */
     //Lsol(ilusch->L, work, y); 
   //} else 
     vbLsolp(nn, ilusch->lu, y, y);//Lsol(ilusch->lu, y, y); 
/*---------------end of SchLsol---------------------------------------
----------------------------------------------------------------------*/
}


void VBSchUsol(vbilutptr ilusch, FLOAT *y) 
{
/*---------------------------------------------------------------------
|  Forward solve for Schur complement part = 
|----------------------------------------------------------------------
| on entry:
| ilusch  = the LU matrix as provided from the ILU functions.
| y       = the right-hand-side vector
|
| on return
| y       = solution of LU x = y. [overwritten] 
|---------------------------------------------------------------------*/
/*-------------------- local variables                        */
  int nn = ilusch->bsz[ilusch->n];//j, *perm = ilusch->rperm,
  //FLOAT *work = ilusch->wk; 
/*-------------------- ONE SIDED ROW PERMS */
   //if (perm != NULL) { 
     //for (j=0; j<n; j++)
       //work[perm[j]] = y[j]; 
/*--------------------  L solve proper */
     //Lsol(ilusch->L, work, y); 
   //} else 
     vbUsolp(nn, ilusch->lu, y, y);//Lsol(ilusch->lu, y, y); 
/*-------------------- case when diagonal scaling is done on columns    */ 
    if (ilusch->D2 !=NULL) 
      vbdscale(nn, ilusch->D2, y, y);

    if (ilusch->dd2 !=NULL) 
      vbdscale(nn, ilusch->dd2, y, y);
/*---------------end of SchLsol---------------------------------------
----------------------------------------------------------------------*/
}


vbp4ptr BLvsol2(FLOAT *x, int nlev, vbp4ptr levmat, vbilutptr ilusch, int flag) 
{
  /* Macro L-solve -- corresponds to left (L) part of arms
  |  preconditioning operation -- 
  |  on entry : 
  |   x =  right- hand side to be operated on by the preconditioner
  |  on return : x is overwritten
  |   x =  output result of operation 
  |  
  |  Note : in-place operation -- b and x can occupy the same space..
  | --------------------------------------------------------------------*/ 
/*-------------------- local variables  */
  int nloc=levmat->bsz[levmat->n], first, lenB;
  vbp4ptr last=levmat; 
/*-------------------- take care of  special cases :  nlev==0 --> lusol  */
  if (nlev == 0) {

     //~ //if (ilusch->lu) 
       VBSchLsol(ilusch, x);//ierr = vblusolC(work,wk,lu);
     //~ //else
	//~ //pSchLUsol(ilusch, x);
       //{SchLsol(ilusch->plu,x);SchUsol(ilusch->plu,x);}

    //VBSchLUsol(ilusch, x);//SchLsol(ilusch,x);
    return (last);
  }
  first = 0;
/*-------------------- descend                                      */
  while (levmat) { 
    nloc = levmat->bsz[levmat->n];
    //printf("\n nloc=%d",nloc); getchar();	
    lenB = levmat->bsz[levmat->nB];
/*-------------------- left scaling                                  */
	//printf("\n levmat->D1=%d,first=%d",levmat->D1,first); getchar();
	//output_dblvector("rhs1.coo" ,&x[first],0, levmat->bsz[levmat->n]);
    if (levmat->D1 !=  NULL) 
      vbdscale(nloc,levmat->D1, &x[first],  &x[first]);//void vbdscale(int nn, double *dd, FLOAT *x, FLOAT * y); 
/*--------------------  RESTRICTION/ DESCENT OPERATION  */
	//printf("\n first=%d",first); getchar();
	//output_dblvector("rhs2.coo" ,&x[first],0, levmat->bsz[levmat->n]);
    if (lenB)
      //output_dblvector("worksecondpart1.coo" ,&x[first],0, 806); 
      vbdescend (levmat, &x[first], &x[first]);//two same input are to avoid the permutation change
    first += lenB; 
    last = levmat;
    levmat = levmat->next;
/*---------------------------------------------------------------------
| next level 
+--------------------------------------------------------------------*/
   }
   if (flag)
   {
  //printf("\n first=%d",first); getchar();
  //~ //if (ilusch->lu) 
     VBSchLsol(ilusch, &x[first]);//ierr = vblusolC(work,wk,lu);
  //~ //else
     //~ //pSchLUsol(ilusch, &x[first]);
     //{SchLsol(ilusch->plu,&x[first]);SchUsol(ilusch->plu,&x[first]);}
  //VBSchLUsol(ilusch, &x[first]);//void VBSchLUsol(vbilutptr ilusch, FLOAT *y)//SchLsol(ilusch,&x[first]);
    }
  return last; 
}

int BUvsol2(FLOAT *x, int nlev, int n, vbp4ptr levmat, vbilutptr ilusch) 
{
  /* Macro U-solve -- corresponds to right (U) part of Barms
  |  preconditioning operation -- 
  |  on entry : 
  |  b  =  right- hand side to be operated on by the preconditioner
  |  on return  = x has been overwritten =
  |  x  =  output result of operation 
  |  
  |  Note : in-place operation -- b and x  can occupy the same space..
  | --------------------------------------------------------------------*/ 

/*-------------------- local variables  */
  int nloc, lenB, first; 
  /*-------------------- work array                                        */
  /*-------------------- take care of  special cases :  nlev==0 --> lusol  */
/*-------------------- case of zero levels                             */
  if (nlev == 0) { 
    VBSchUsol(ilusch, x);  
    return(0);
  }
/*-------------------- general case                               */
  nloc = levmat->bsz[levmat->n]; 
  lenB = levmat->bsz[levmat->nB]; 
  first = n - nloc; 
/*-------------------- last level                                 */
  first += lenB; 
//printf("\n first=%d",first); getchar();
  VBSchUsol(ilusch, &x[first]); 
/*-------------------- other levels                               */
  while (levmat) {
    nloc = levmat->bsz[levmat->n]; 
    first -= levmat->bsz[levmat->nB];
	//printf("\n first=%d",first); getchar();
    if (levmat->n) 
      vbascend(levmat, &x[first],&x[first]);
/*-------------------- right scaling */
    //printf("\n levmat->D2=%d,first=%d",levmat->D2,first); getchar();
    if (levmat->D2 !=  NULL) 
      vbdscale(nloc, levmat->D2, &x[first], &x[first]) ;//void vbdscale(int n, double *dd, FLOAT *x, FLOAT * y, int *bsz)
    levmat = levmat->prev; 
  }
  return 0;
/*--------------------  PROLONGATION/ ASCENT OPERATION */
}

int barmsol2(FLOAT *x,  vbarms Prec) 
{ /* combined preconditioning operation -- combines the
  |  left and right actions. 
  | 
  |  on entry : 
  |   x =  right- hand side to be operated on by the preconditioner
  |  on return : x is overwritten - 
  |   x =  output result of operation 
  |  
  |  Note : in-place operation -- b and x can occupy the same space..
  | --------------------------------------------------------------------*/ 
/*-------------------- local variables  */
  vbp4ptr levmat = Prec->levmat;
  vbilutptr ilusch = Prec->ilus;
  int nlev = Prec->nlev;
  int n; 
  vbp4ptr last;
  if (nlev == 0) {
    n = ilusch->bsz[ilusch->n];
     if (ilusch->lu) 
       VBSchLUsol(ilusch, x);//ierr = vblusolC(work,wk,lu);
     else
       pSchLUsol(ilusch, x);//void pSchLUsol(vbilutptr ilusch, FLOAT *y) //{SchLsol(ilusch->plu,x);SchUsol(ilusch->plu,x);}
    //VBSchLUsol(ilusch, x);
    //SchLsol(ilusch, x);//pliut 
    //SchUsol(ilusch, x); 
    return 0;
  }
  else
    n = levmat->bsz[levmat->n];
  last = BLvsol2(x, nlev, levmat, ilusch, 1) ;
  BUvsol2(x, nlev, n, last, ilusch) ; 
  return 0; 
}
