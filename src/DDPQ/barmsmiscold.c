#include "protos.h"

void swapmm(BData v[], int i, int j) {//for dense columns switch
    BData temp;
    temp = v[i];
    v[i] = v[j];
    v[j] = temp;
}

double vbnorm2( int sz, FLOAT *a );

int vbplussvb1(int job, vbsptr A, FLOAT s, vbsptr B, vbsptr C, double tol){
    /*  C= A + s B (not sorted), with dropping strategy
NOTE: 
    No require the elements of each row are sorted with increasing column indices in
    each row.

*/

    int i, ii, ka, kb, values, len, jcol, jpos,dimR, dimC, blocksz, inc = 1,len1;
    int *iw, *ja, *jb, *jc, *jcdrop, nrow, ncol,*nB,*bsz = A->bsz;
    BData *ba, *bb, *w, *wdrop;
    double tnorm,t,tolnorm;

    nrow = A->n;////
    ncol = A->n;////
    iw =(int *)Malloc(sizeof(int)*ncol,"vbplussvb iw");//memory n later
    jc =(int *)Malloc(sizeof(int)*ncol,"vbplussvb jc");//memory n later
    jcdrop =(int *)Malloc(sizeof(int)*ncol,"vbplussvb jc");

    w  =(BData *)Malloc(sizeof(BData)*ncol,"vbplussvb work array");
    wdrop =(BData *)Malloc(sizeof(BData)*ncol,"vbplussvb work array");

    //init
    for(i=0; i<ncol; i++)//not beyond the column number of B
        iw[i] = -1;

    nB =(int *)Malloc(sizeof(int)*nrow,"vbplussvb nB");//memory n later
    for(i=0;i<nrow;i++)
    {
        nB[i] = B_DIM(bsz,i);
    }

    if (setupVBMat(C,nrow,nB)) return 1;
    free(nB);

    /*#ifdef DOUBLEELETYPE
    setupCS(C, nrow, 1);
#else
    zsetupCS(C, nrow);
#endif*/

    values = (job!=0);
    //	printf("\n vbplussvb 2"); getchar();
    //LOOP row of A
    for(i=0; i< nrow; i++){
        //printf("\n vbplussvb 3,row=%d",i); getchar();
        dimR = B_DIM(bsz,i);
        len = -1;////
        ba = A->ba[i];//pointer
        ja = A->ja[i];//pointer
        for(ka=0; ka< A->nzcount[i]; ka++){
            len ++;
            jcol = ja[ka];
            /*if (len > nzmax)
                return 1;*/
            jc[len] = jcol;
            if (values) {//w[len] = ba[ka];
                dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
                w[len] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                copyBData( dimR, dimC, w[len], ba[ka], 0 );
                //output_dblvector("w[len].coo",w[len],0, blocksz);getchar();
            }
            iw[jcol] = len;
        }//for ka=0
        bb = B->ba[i];
        jb = B->ja[i];
        for(kb=0; kb<B->nzcount[i]; kb++){
            jcol = jb[kb];
            jpos = iw[jcol];
            dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
            if (jpos < 0){//not jpos ==0
                len++;
                /*if(len >nzbax)
                    return 2;*/
                jc[len] = jcol;
                if (values) {//w[len] = s*bb[kb];
                    //dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
                    //printf("\n in 3, dimC=%d,dimR=%d",dimC,dimR); getchar();
                    w[len] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                    //printf("\n in 3, w[kc]=%d,bb[kb][i]=%f",w[kc],bb[kb][i]); getchar();
                    for( ii = 0; ii < blocksz; ii++ ) w[len][ii] = s*bb[kb][ii];//copyBData( dimR, dimC, w[kc], ba[ka], 1 );
                    //output_dblvector("ba[ka].coo" ,ba[ka],0,blocksz );
                }
                iw[jcol] = len;
            }
            else{
                if(values) {//w[jpos] += s*bb[kb];
                    GAXPY(blocksz,s,bb[kb],inc,w[jpos],inc);//w[kc]= ba[ka] + s*bb[kb];work well
                }
            }//if jpos==0
        }//for kb=0
        //copy to matrix C
        tnorm = 0;
        for (ii=0; ii<=len; ii++){//dropping strategy
            jcol = jc[ii];
            dimC = B_DIM(bsz,jcol);
            t = vbnorm2( dimR*dimC, w[ii] );
            tnorm = max(t,tnorm);
        }
        len1 = 0;tolnorm = tol*tnorm;
        for (ii=0; ii<=len; ii++){//dropping strategy
            jcol = jc[ii];
            dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
            if (vbnorm2(blocksz,w[ii]) > tolnorm){//vbnorm2(blocksz,w[ii]) > tol//vbnorm2( dim * szjrow, buf_fact ) * xnrm[jrow] <= tolnorm//w[j] = DNRM2(blocksz,ba[j],inc);
                wdrop[len1] = w[ii];//copyBData( dimR, dimC, wdrop[len1], w[ii], 0 );//wdrop[len1] = w[ii];void copyBData( int m, int n, BData dst, BData src, int isig )
                jcdrop[len1] = jc[ii];//jcdrop[len1] = jc[ii];
                len1++;//len1++;
            }
            else free(w[ii]);
        }
        /*if (len>=0){
            len ++;
            C->ja[i] = (int *)Malloc(sizeof(int)*len,"vbplussvb C->ja");
            C->ba[i] = (BData *)Malloc(sizeof(BData)*len,"vbplussvb c->ba");
            memcpy(C->ja[i], jc, len*sizeof(int));
            memcpy(C->ba[i],  w, len*sizeof(BData));
            C->nzcount[i] = len;
        }
        else{
            C->ja[i] = NULL;
            C->ba[i] = NULL;
            C->nzcount[i] = 0;
        }*/
        //printf("\n len=%d,len1=%d",len,len1); getchar();
        if (len1>0){
            len ++;
            C->ja[i] = (int *)Malloc(sizeof(int)*len1,"vbplussvb C->ja");//????
            C->ba[i] = (BData *)Malloc(sizeof(BData)*len1,"vbplussvb c->ba");
            memcpy(C->ja[i], jcdrop, len1*sizeof(int));
            memcpy(C->ba[i],  wdrop, len1*sizeof(BData));
            C->nzcount[i] = len1;
        }
        else{
            C->ja[i] = NULL;
            C->ba[i] = NULL;
            C->nzcount[i] = 0;
        }
        for (ii=0; ii< len; ii++)
            iw[jc[ii]] = -1;////
    }//for i=0
    //	printf("\n vbplussvb 4"); getchar();
    free(w);
    free(wdrop);
    free(iw);
    free(jc);
    free(jcdrop);
    //	printf("\n vbplussvb 5"); getchar();
    return 0;
}

/*int vbplussvb1(int job, vbsptr A, double s, vbsptr B, vbsptr C){

    int i, ii, ka, kb, values, nzmax, len, jcol, jpos,dimR, dimC, blocksz, inc = 1,len1;
    int *iw, *ja, *jb, *jc,*jcdrop, nrow, ncol,*nB,*bsz = A->bsz;
    double tol = 0.00001;
    BData *ba, *bb, *w, *wdrop;

    nrow = A->n;////
    ncol = A->n;////
    iw =(int *)Malloc(sizeof(int)*ncol,"vbplussvb iw");//memory n later
    jc =(int *)Malloc(sizeof(int)*ncol,"vbplussvb jc");//memory n later
    jcdrop =(int *)Malloc(sizeof(int)*ncol,"vbplussvb jc");//memory n later

    w  =(BData *)Malloc(sizeof(BData)*ncol,"vbplussvb work array");
    wdrop  =(BData *)Malloc(sizeof(BData)*ncol,"vbplussvb work array");
    nzmax = nrow;////default
    //init
    for(i=0; i<ncol; i++)//not beyond the column number of B
        iw[i] = -1;

    nB =(int *)Malloc(sizeof(int)*nrow,"vbplussvb nB");//memory n later
        for(i=0;i<nrow;i++)
        {
            nB[i] = B_DIM(bsz,i);
        }

        if (setupVBMat(C,nrow,nB)) return 1;
    free(nB);


    values = (job!=0);
//	printf("\n vbplussvb 2"); getchar();
    //LOOP row of A
    for(i=0; i< nrow; i++){
    //printf("\n vbplussvb 3,row=%d",i); getchar();
        dimR = B_DIM(bsz,i);
        len = -1;////
        ba = A->ba[i];//pointer
        ja = A->ja[i];//pointer
        for(ka=0; ka< A->nzcount[i]; ka++){
            len ++;
            jcol = ja[ka];

            jc[len] = jcol;
            if (values) {//w[len] = ba[ka];
            dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
            w[len] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
            copyBData( dimR, dimC, w[len], ba[ka], 0 );
            //output_dblvector("w[len].coo",w[len],0, blocksz);getchar();
            }
            iw[jcol] = len;
        }//for ka=0
        bb = B->ba[i];
        jb = B->ja[i];
        for(kb=0; kb<B->nzcount[i]; kb++){
            jcol = jb[kb];
            jpos = iw[jcol];
            dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
            if (jpos < 0){//not jpos ==0
                len++;
                jc[len] = jcol;
                if (values) {//w[len] = s*bb[kb];
                      //dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
                      //printf("\n in 3, dimC=%d,dimR=%d",dimC,dimR); getchar();
                      w[len] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                      //printf("\n in 3, w[kc]=%d,bb[kb][i]=%f",w[kc],bb[kb][i]); getchar();
                      for( ii = 0; ii < blocksz; ii++ ) w[len][ii] = s*bb[kb][ii];//copyBData( dimR, dimC, w[kc], ba[ka], 1 );
                      //output_dblvector("ba[ka].coo" ,ba[ka],0,blocksz );
                }
                iw[jcol] = len;
            }
            else{
                if(values) {//w[jpos] += s*bb[kb];
                      GAXPY(blocksz,s,bb[kb],inc,w[jpos],inc);//w[kc]= ba[ka] + s*bb[kb];work well
                }
            }//if jpos==0
        }//for kb=0
        //copy to batrix C
        len1 = 0;
        for (ii=0; ii<=len; ii++){//dropping strategy
            jcol = jc[ii];
            dimC = B_DIM(bsz,jcol);blocksz = dimR*dimC;
            if (vbnorm2(blocksz,w[ii]) > tol){//vbnorm2(blocksz,w[ii]) > tol//vbnorm2( dim * szjrow, buf_fact ) * xnrm[jrow] <= tolnorm//w[j] = DNRM2(blocksz,ba[j],inc);
                wdrop[len1] = w[ii];//copyBData( dimR, dimC, wdrop[len1], w[ii], 0 );//wdrop[len1] = w[ii];void copyBData( int m, int n, BData dst, BData src, int isig )
                jcdrop[len1] = jc[ii];//jcdrop[len1] = jc[ii];
                len1++;//len1++;
                }
//			else free(w[ii]);
        }
        if (len1>0){
            C->ja[i] = (int *)Malloc(sizeof(int)*len1,"vbplussvb C->ja");
            C->ba[i] = (BData *)Malloc(sizeof(BData)*len1,"vbplussvb c->ba");
            memcpy(C->ja[i],jcdrop, len1*sizeof(int));
            memcpy(C->ba[i], wdrop, len1*sizeof(BData));
            C->nzcount[i] = len1;
        }
        else{
            C->ja[i] = NULL;
            C->ba[i] = NULL;
            C->nzcount[i] = 0;
        }
        for (ii=0; ii< len; ii++)
            iw[jc[ii]] = -1;////
    }//for i=0
//	printf("\n vbplussvb 4"); getchar();
    free(w);
    free(iw);
    free(jc);
    free(jcdrop);
    free(wdrop);
//	printf("\n vbplussvb 5"); getchar();
    return 0;
}*/



/*int vbplussvb1(int job, vbsptr A, FLOAT s, vbsptr B, vbsptr C){
//  C= A + s B (sorted)(vbsptr)
NOTE: 
    Require the elements of each row are sorted with increasing column indices in
    each row.


    int i,j1, j2, ka, kb, kc, values, nzmax, nrow, ncol, dimR, dimC, blocksz, inc = 1, ii;
    int *iw, *ja, *jb, *jc, kamax, kbmax,*nB,*bsz = A->bsz;
    BData *ba,*bb,*w;

    nrow = A->n;
    ncol = A->n;
    //printf("\n vbplussvb1, An=%d,Bn=%d", A->n,B->n); getchar();
    iw =(int *)Malloc(sizeof(int)*nrow,"vbplussvb1 iw");//memory n later
    jc =(int *)Malloc(sizeof(int)*ncol,"vbplussvb1 jc");//memory n later

    w    = (BData *)Malloc(sizeof(BData)*nrow,"vbplussvb1 work array");
    nzmax = nrow;////default

    //init
    for(i=0; i< nrow; i++)//not beyond the column number of B
        iw[i] = -1;

    nB =(int *)Malloc(sizeof(int)*nrow,"vbplussvb1 nB");//memory n later
        for(i=0;i<nrow;i++)
        {
            nB[i] = B_DIM(bsz,i);
        }

        if (setupVBMat(C,nrow,nB)) return 1;
    free(nB);

    values = (job!=0);

    //LOOP row of A
    for(i=0; i< nrow; i++){
        ka = 0;
        kb = 0;
        kc = 0;
        kamax = A->nzcount[i];
        kbmax = B->nzcount[i];
        ba = A->ba[i];
        bb = B->ba[i];
        ja = A->ja[i];
        jb = B->ja[i];
        dimR = B_DIM(bsz,i);
        //printf("\n vbplussvb1 i=%d",i); getchar();
        //printf("\n vbplussvb1 kamax=%d,kbmax=%d",kamax,kbmax); getchar();
        while((ka<kamax) || (kb<kbmax)){
            if (ka < kamax)
                j1 = ja[ka];
            else
                j1 = ncol + 1;
            if (kb < kbmax)
                j2 = jb[kb];
            else
                j2 = ncol + 1;
            //three cases
            //printf("\n vbplussvb1 j1=%d,j2=%d",j1,j2); getchar();
            //printf("\n before  ka=%d,kb=%d,kc=%d",ka,kb,kc); getchar();
            if (j1 == j2){
                if (values) {
                      dimC = B_DIM(bsz,j1);blocksz = dimR*dimC;
                      w[kc] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                      copyBData( dimR, dimC, w[kc], ba[ka], 0 );
                      GAXPY(blocksz,s,bb[kb],inc,w[kc],inc);//w[kc]= ba[ka] + s*bb[kb];work well
                }
                jc[kc] = j1;
                ka ++;
                kb ++;
                kc ++;
            }
            else if (j1 < j2){
                jc[kc] = j1;
    //printf("\n vbplussvb1 ka=%d,kb=%d,kc=%d",ka,kb,kc); getchar();
                if (values){ //w[kc] = ba[ka];
                      dimC = B_DIM(bsz,j1);blocksz = dimR*dimC;
                      w[kc] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                      copyBData( dimR, dimC, w[kc], ba[ka], 0 );
                }
                ka ++;
                kc ++;
            }
            else if (j1 > j2){//j1 > j2
                jc[kc] = j2;
    //printf("\n in 3, ka=%d,kb=%d,kc=%d",ka,kb,kc); getchar();
                if (values){ //w[kc] = s*bb[kb];
                      dimC = B_DIM(bsz,j2);blocksz = dimR*dimC;
                      //printf("\n in 3, dimC=%d,dimR=%d",dimC,dimR); getchar();
                      w[kc] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                      //printf("\n in 3, w[kc]=%d,bb[kb][ii]=%f",w[kc],bb[kb][ii]); getchar();
                      for( ii = 0; ii < blocksz; ii++ ) w[kc][ii] = s*bb[kb][ii];//copyBData( dimR, dimC, w[kc], ba[ka], 1 );
                      //output_dblvector("ba[ka].coo" ,ba[ka],0,blocksz );
                }
                kb ++;
                kc ++;
            }
            //if( kc> nzmax) return 1;
        //printf("\n vbplussvb1 ka=%d,kb=%d,kc=%d",ka,kb,kc); getchar();
        }//while
        //copy to matrix C

        if (kc > 0){
            //no need to do kc++;
            C->ja[i] = (int *)Malloc(sizeof(int)*(kc),"vbplussvb1 C->ja");
            C->ba[i] = (BData *)Malloc(sizeof(BData)*(kc),"vbplussvb1 c->ma");
            memcpy(C->ja[i], jc, (kc)*sizeof(int));
            memcpy(C->ba[i],  w, (kc)*sizeof(BData));
            C->nzcount[i] = kc;
        }
        else{
            C->ja[i] = NULL;
            C->ba[i] = NULL;
            C->nzcount[i] = 0;
        }
    }//for ii=0

    free(w);
    free(iw);
    free(jc);
    return 0;
} */

int vbmulvb(int job, vbsptr A, vbsptr B, vbsptr C){
    /*
   C= A*B

*/
    int i, ka, kb, ii, jj, values, len, jcol, jpos, dimR, dimC, BdimC, blocksz, inc = 1;
    int *iw,*ja,*jb,*jc, nrow, ncol,*bsz = A->bsz, *bszc = A->bszc, *nB;
    int max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(FLOAT);
    FLOAT one = 1.0, zero = 0.0;
    BData *ba,*bb,*w = NULL, scal = NULL,workbuff=NULL;
    nrow = A->n;
    ncol = max(A->n,B->n);//just for allocating memory(not means the number of columns of some matrix)
    //printf("\n vbmulvb max_blk_sz=%d",max_blk_sz); getchar();
    workbuff = (BData)Malloc( max_blk_sz, "vbmulvb" );//workbuff = (BData)Malloc( max_blk_sz, "vbmulvb" );

    //nzmax = n;////??
    //printf("\n vbmulvb begin nrow=%d, ncol=%d",nrow,ncol); getchar();

    if (ncol<=0) {printf("error in vbmulvb"); getchar();}

    iw =(int *)Malloc(sizeof(int)*ncol,"vbmulvb iw");//memory n later

    //printf("\n vbmulvb begin 22-1"); getchar();
    jc =(int *)Malloc(sizeof(int)*ncol,"vbmulvb jc");//memory n later

    //printf("\n vbmulvb begin 22-2"); getchar();
    w  =(BData *)Malloc(sizeof(BData)*ncol,"vbmulvb work array");


    //	printf("\n vbmulvb begin 22"); getchar();
    //init
    for(i=0; i<A->n; i++)//not beyond the column number of B
        iw[i] = -1;

    nB =(int *)Malloc(sizeof(int)*nrow,"vbmulvb nB");//memory n later
    for(i=0;i<nrow;i++)
    {
        nB[i] = B_DIM(bsz,i);
    }

    if (setupVBMat(C,nrow,nB)) return 1;
    free(nB);nB=NULL;

    /*#ifdef DOUBLEELETYPE
    setupCS(C, nrow, 1);//
#else
    zsetupCS(C, nrow);
#endif*/
    //	printf("\n vbmulvb begin 33"); getchar();
    values = (job!=0);

    //LOOP row of A
    for(ii=0; ii< nrow; ii++){
        dimR = B_DIM(bsz,ii);
        ba = A->ba[ii];//pointer
        ja = A->ja[ii];//pointer
        len = -1;////

        //printf("\n vbmulvb nnz(A[%d])=%d",ii,A->nzcount[ii]);getchar();
        for(ka=0; ka< A->nzcount[ii]; ka++){
            dimC = B_DIM(bszc,ja[ka]);
            if (values) scal = ba[ka];//??
            jj = ja[ka];

            bb = B->ba[jj];//pointer 2009.07.29 move here
            jb = B->ja[jj];//pointer
            //printf("\n vbmulvb nnz(B[%d])=%d",jj,B->nzcount[jj]);getchar();
            for(kb=0; kb<B->nzcount[jj]; kb++){
                //mb = B->ma[jj];//pointer
                //jb = B->ja[jj];//pointer
                BdimC = B_DIM(bsz,jb[kb]);//here, bsz is also the block column index of B,
                jcol = jb[kb];
                jpos = iw[jcol];
                blocksz = dimR*BdimC;
                //printf("\n vbmulvb jcol=%d",jcol);getchar();
                if (jpos == -1){////already value copied or not,the colum index of C, iw is the array to record this.
                    len = len+1;
                    /*if (len > nzmax){//nzmax later
                        return 1;
                    }*/
                    jc[len] = jcol;//store column index of C
                    iw[jcol] = len;

                    if (values){
                        w[len] = (BData)Malloc( blocksz*sizeof(FLOAT), "vbmulvb" );//w[len] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));		//w[len] = scal*bb[kb];//??
                        //w[len] = (FLOAT*)realloc(w[len],blocksz*sizeof(FLOAT));
                        //printf("\n dimR=%d,dimC=%d,BdimC=%d,len=%d",dimR,dimC,BdimC,len); getchar();
                        //output_dblvector("bb[kb].coo" ,bb[kb],0, 16);getchar();
                        //output_dblvector("scal.coo" ,scal,0, 16);getchar();
                        GGEMM ("n","n",dimR, BdimC, dimC, one, scal, dimR, bb[kb], dimC, zero, w[len], dimR );//GGEMM ("n","n",  dim, sz, dim, one,D[i], dim, ba[j], dim, zero, buf,dim) ;
                        //output_dblvector("w[len].coo",w[len],0, 16);getchar();
                    }
                }
                else{
                    if (values){
                        //printf("\n dimR=%d,dimC=%d,BdimC=%d,jpos=%d,in +",dimR,dimC,BdimC,jpos); getchar();
                        //output_dblvector("bb[kb]+.coo" ,bb[kb],0, 16);getchar();
                        //output_dblvector("scal+.coo" ,scal,0, 16);getchar();
                        //workbuff = (FLOAT*)realloc(workbuff,blocksz*sizeof(FLOAT));//workbuff = (FLOAT*)malloc(blocksz*sizeof(FLOAT));//nB =(int *)Malloc(sizeof(int)*nrow,"vbmulvb nB");
                        // workbuff = (BData)Malloc( blocksz*sizeof(FLOAT), "vbmulvb" );
                        GGEMM ("n","n",dimR, BdimC, dimC, one, scal, dimR, bb[kb], dimC, zero, workbuff, dimR );
                        //output_dblvector("workbuff+.coo" ,workbuff,0, 16);getchar();
                        //output_dblvector("w[jpos]pre+.coo" ,w[jpos],0, 16);getchar();
                        GAXPY(blocksz,one,workbuff,inc,w[jpos],inc);//for(kk=0;kk<blocksz;kk++) w[jpos][kk] += workbuff[kk];//GAXPY(blocksz,one,workbuff,inc,w[jpos],inc);
                        //free(workbuff);
                        //output_dblvector("w[jpos]+.coo",w[jpos],0, 16);getchar();
                    }						//w[jpos] += scal*bb[kb];//??//??
                }//if jpos==0
            }//for kb
            //free??
        }//for ka=0
        //copy to matrix C
        if (len>=0){
            len++;//important
            C->ja[ii] = (int *)Malloc(sizeof(int)*len,"vbmulvb C->ja");
            C->ba[ii] =(BData *)Malloc(sizeof(BData)*len,"vbmulvb C->ba");
            memcpy(C->ja[ii], jc, len*sizeof(int));
            memcpy(C->ba[ii],  w, len*sizeof(BData));
            C->nzcount[ii] = len;
        }
        else{
            C->ja[ii] = NULL;
            C->ba[ii] = NULL;
            C->nzcount[ii] = 0;
        }
        //output_intvector("iw.coo" ,iw,0, ncol); getchar();
        for (i=0; i<len; i++)
            iw[jc[i]] = -1;//?
    }//for ii=0

    //printf("\nwarning: no free in vbmulvb");getchar();
    //outputvbmat1(C,"Cintimes.coo",1);
    free(w);	free(iw);	free(jc);    free(workbuff);
    //printf("\nwarning: no free in vbmulvb");getchar();
    return 0;
}


int setupVBvector( vbbsptr vbbmat, int nset, int *nBB)
{
    /*----------------------------------------------------------------------
| Initialize VBBSpaFmt structs
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbbmat ) =  Pointer to a VBBSpaFmt struct.
|      
|
| On return:
|===========
|
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    /*typedef struct VBBSpaFmt {
    int n;
    int *bsz;
    vbsptr *vba;
} VBBSparMat, *vbbsptr;*/
    int i;
    vbbmat->n = nset;
    //vbbmat->nc = nc;
    //printf("\n vbn=%d,vbnc=%d",vbbmat->n,vbbmat->nc); getchar();
    if( nBB ) {
        vbbmat->bsz = (int *)Malloc( sizeof(int)*(nset+1), "setupvbbmat" );
        vbbmat->bsz[0] = 0;
        for( i = 1; i <= nset; i++ ) {
            vbbmat->bsz[i] = vbbmat->bsz[i-1] + nBB[i-1];
            //printf("\n vbn=%d",vbbmat->bsz[i]); getchar();
        }
    } else
        vbbmat->bsz = NULL;
    
    vbbmat->vba = (vbsptr*)Malloc( sizeof(vbsptr) * nset, "setupvbbmat" );
    //printf("\n vbn=%d,vbnc=%d",vbbmat->n,vbbmat->nc); getchar();
    return 0;

}


int setupVBilutvector( vbbiluptr vbbmat, int nset, int *nBB)
{
    /*----------------------------------------------------------------------
| Initialize VBBSpaFmt structs
|----------------------------------------------------------------------
| on entry:
|==========
| ( vbbmat ) =  Pointer to a VBBSpaFmt struct.
|      
|
| On return:
|===========
|
|
| integer value returned:
|             0   --> successful return.
|            -1   --> memory allocation error.
|--------------------------------------------------------------------*/
    /*typedef struct VBBILUfac {
    int n;
    int *bsz;
    vbiluptr *vlu;
} VBBILUSpar, *vbbiluptr;*/
    int i;
    vbbmat->n = nset;
    //vbbmat->nc = nc;
    //printf("\n vbn=%d,vbnc=%d",vbbmat->n,vbbmat->nc); getchar();
    if( nBB ) {
        vbbmat->bsz = (int *)Malloc( sizeof(int)*(nset+1), "setupvbbmat" );
        vbbmat->bsz[0] = 0;
        for( i = 1; i <= nset; i++ ) {
            vbbmat->bsz[i] = vbbmat->bsz[i-1] + nBB[i-1];
            //printf("\n vbn=%d",vbbmat->bsz[i]); getchar();
        }
    } else
        vbbmat->bsz = NULL;
    
    vbbmat->vlu = (vbiluptr*)Malloc( sizeof(vbiluptr) * nset, "setupvbbmat" );
    //printf("\n vbn=%d,vbnc=%d",vbbmat->n,vbbmat->nc); getchar();
    return 0;

}




int Bsplit(vbsptr B, vbbsptr Bb, vbsptr F, vbbsptr Ff, int *nBB, int nset)
{

    /*typedef struct VBBSpaFmt {
    int n;
    int *bsz;
    vbsptr *vba;
} VBBSparMat, *vbbsptr;*/
    /*typedef struct VBSpaFmt {
    int n;
    int *bsz;
    int *nzcount;
    int **ja;
    BData **ba;
    BData *D;
} VBSparMat, *vbsptr;*/     

    int j, i,k,*ja,*jaF, *nB, *nC, dimR, dimC,firstBi,rowBi,nzcount,nzcountF, col,*Bija,*Fija,maxsize = 0;
    int *bsz = B->bsz,*bszz,*bszc = F->bszc;
    BData *ba,*Biba,*baF,*Fiba;
    vbsptr *vba,*vbaF;
    vbsptr Bi,Fi;


    // Biba = (BData *) Malloc(bsize*sizeof(BData), "vbSplit4copy:3" );

    if (setupVBvector(Bb,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
    //nBB is size of diagonal coarse block, nset is the dimension of it.
    if (setupVBvector(Ff,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)

    //Bb->n = nset;Ff->n = nset;
    //printf("\nafter setup "); getchar();
    bszz = Bb->bsz;//bszz is the coarse block index

    nB = (int *) Malloc(B->n*sizeof(int), "vbSplit4copy:0" );
    nC = (int *) Malloc(F->nc*sizeof(int), "vbSplit4copy:0" );
    for (i=0; i<B->n; i++){
        nB[i]=B_DIM(bsz,i);
    }//nB is size of diagonal fine block, B->n is the dimension of it.
    //output_intvector("nBinB.coo" ,nB,0, B->n);getchar();
    for (i=0; i<F->nc; i++){
        nC[i]=B_DIM(bszc,i);
    }
    //printf("\nafter setup "); getchar();
    for( i = 0; i < nset; i++ ) {
        //printf("\ni=%d",i); getchar();
        vba = Bb->vba;
        vbaF = Ff->vba;
        Bi = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
        Fi = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
        firstBi = bszz[i];//in the fine block vector, find the first address of each coarse block
        //printf("\nfirstBi=%d",firstBi); getchar();
        if (setupVBMat(Bi,nBB[i],&nB[firstBi])) return 1; //if (setupVBMat(B,bsize,nB)) goto label111;
        //if (setupVBMat1(F,bsize,csize,nB,nC)) goto label111;
        if (setupVBMat1(Fi,nBB[i],F->nc,&nB[firstBi],nC)) return 1; //int setupVBMat1( vbsptr vbmat, int n,int nc, int *nBR, int *nBC )
        maxsize = max(maxsize,Bi->bsz[Bi->n]);
        //output_intvector("bszinBi.coo" ,Bi->bsz,0, nBB[i]+1);getchar();
        //printf("\Bi->n=%d",Bi->n); getchar();

        for( j = 0; j < nBB[i]; j++ ) {
            rowBi = firstBi+j;
            dimR = B_DIM(bsz,rowBi);
            nzcount = Bi->nzcount[j] = B->nzcount[rowBi];
            nzcountF = Fi->nzcount[j] = F->nzcount[rowBi];
            ja = B->ja[rowBi];
            ba = B->ba[rowBi];
            jaF = F->ja[rowBi];
            baF = F->ba[rowBi];
            if (nzcount>0) {
                Bi->ja[j] = (int *) Malloc(nzcount*sizeof(int), "vbSplit4copy:5" );
                Bi->ba[j] = (BData *) Malloc(nzcount*sizeof(BData), "vbSplit4copy:6" );
            }
            if (nzcountF>0) {
                Fi->ja[j] = (int *) Malloc(nzcountF*sizeof(int), "vbSplit4copy:5" );
                Fi->ba[j] = (BData *) Malloc(nzcountF*sizeof(BData), "vbSplit4copy:6" );
            }
            //printf("\ni=%d,j=%d,nz=%d",i,j,nzcount); getchar();
            // memcpy(Bi->ja[j], ja, nzcount*sizeof(int));//ja-
            //printf("\ni=%d,j=%d",i,j); getchar();
            Bija = Bi->ja[j];
            Biba = Bi->ba[j];
            Fija = Fi->ja[j];
            Fiba = Fi->ba[j];
            for( k = 0; k < nzcount; k++ ) {
                //printf("\ni=%d,j=%d,rowBi=%d,k=%d,col=%d",i,j,rowBi,k,col); getchar();
                col = ja[k];
                dimC = B_DIM(bsz,col);
                Bija[k] = ja[k]- firstBi;
                Biba[k] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
                copyBData( dimC, dimR, Biba[k], ba[k], 0 );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
            }
            for( k = 0; k < nzcountF; k++ ) {
                //printf("\ni=%d,j=%d,rowBi=%d,k=%d,col=%d",i,j,rowBi,k,col); getchar();
                col = jaF[k];
                dimC = B_DIM(bszc,col);
                Fija[k] = jaF[k];
                Fiba[k] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
                copyBData( dimC, dimR, Fiba[k], baF[k], 0 );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
            }
        }
        vba[i] = Bi;
        vbaF[i] = Fi;
        //printf("\n vba[i]->n=%d",vba[i]->n); getchar();
    }
    Bb->maxsize = maxsize;
    Ff->maxsize = maxsize;
    // printf("\n Bb->maxsize=%d",Bb->maxsize); getchar();
    if (nB) free(nB);
    if (nC) free(nC);
    return 0;
}



void qqsortm(int *ja,int *jab, BData *ma, int left, int right){//new cs2dptr
    /*----------------------------------------------------------------------
|
| qqsort: sort ja[left]...ja[right], and ma[left]...ma[right] into increasing order
| from Kernighan & Ritchie
|
| ma holds the real values
|
|---------------------------------------------------------------------*/
    int i, last;
    if (left >= right)  return;
    swapj(ja, left, (left+right)/2);
    swapj(jab, left, (left+right)/2);
    swapmm(ma, left, (left+right)/2);
    last = left;
    for (i=left+1; i<=right; i++) {
        if (ja[i] < ja[left]) {
            swapj(ja, ++last, i);
            swapj(jab, last, i);
            swapmm(ma, last, i);
        }
    }
    swapj(ja, left, last);
    swapj(jab, left, last);
    swapmm(ma, left, last);
    qqsortm(ja, jab, ma, left, last-1);
    qqsortm(ja, jab, ma, last+1, right);
}

int vb2csr2d(vbsptr vbmat, cs2dptr mat)
{
    /*
 * This  subroutine converts the vbsptr matrix to sparse 2D FLOAT pointer format
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 *
 * vbmat = Various Block Sparse Row format Matrix, vbsptr

typedef struct Spa2Dmt {
--------------------------------------------- 
| sparse 2D pointer to matrix,
| only store the nonzero columns. 
| specially for F in Barms
|---------------------------------------------
  int n;
  int nc;
  int nzcount;
  int *ja;
  int *jab
  FLOAT **ma;
} Spar2DMat, *cs2dptr;
 *
 * on return:
 *-----------
 *
 *  2D sparse FLOAT pointer format matrix
 *
 *  ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int i, j, k, nzcount, col,  dimR, dimC, len = -1, firstrow, jpos, colp;
    int n = vbmat->n,nc=vbmat->nc, *ja, *bsz = vbmat->bsz,*bszc=vbmat->bszc,nnr=bsz[n],nnc=bszc[nc],*iw,*jc,*jcb;
    //FLOAT *w;
    BData *ba,*w;

    mat->n = nnr;
    mat->nc = nnc;
    mat->nzcount = 0;
    mat->ja = NULL;
    mat->jab = NULL;
    mat->ma = NULL;
    //for(i=0;i<nnc;i++)
    //{
    //mat[i]=(FLOAT*)malloc(nnr*sizeof(FLOAT));
    //memset(mat[i], 0, nnr*sizeof(FLOAT));
    //}
    //setupCS(csmat,n,1);
    //printf("\n n1=%d",csmat->n); getchar();

    //w =(FLOAT*)malloc(n*sizeof(FLOAT));

    w =(FLOAT**)malloc(nnc*sizeof(FLOAT*));//reallocate later
    iw =(int*)malloc(nnc*sizeof(int));
    jc =(int*)malloc(nnc*sizeof(int));
    jcb = (int*)malloc(nnc*sizeof(int));

    for(i=0; i<nnc; i++)//
        iw[i] = -1;

    for( i = 0; i < n; i++ )
    {
        //printf("\n i=%d",i); getchar();
        dimR = B_DIM(bsz,i);//get the row dimension of the block
        nzcount = vbmat->nzcount[i];//
        //printf("\n csmat->nzcount[i]=%d",csmat->nzcount[i]); getchar();
        ja = vbmat->ja[i];
        //printf("\n vbmat->ja[i]=%d",vbmat->ja[i]); getchar();
        ba = vbmat->ba[i];
        firstrow = bsz[i];
        for( j = 0; j < nzcount; j++ )
        {
            //printf("\n j=%d nzcount=%d",j,nzcount); getchar();
            col = ja[j];
            dimC = B_DIM(bszc,col);//obtain the column dimension of the block
            //blocksz=dimR*dimC;//caculate the blocksize
            for(k=0; k< dimC;k++)
            {
                colp = bszc[col]+k;
                //jcol = jb[kb];
                jpos = iw[colp];
                if (jpos == -1){////already value copied or not,the colum index of C, iw is the array to record this.
                    len = len+1;
                    /*if (len > nzmax){//nzmax later
                return 1;
            }*/
                    jc[len] = colp;//store column index of C
                    jcb[len] = col;//store block column index of C
                    iw[colp] = len;
                    w[len] = (BData)Malloc( nnr*sizeof(FLOAT), "vb2csr2d" );//w[len] = (FLOAT*)malloc(blocksz*sizeof(FLOAT));		//w[len] = scal*bb[kb];//??//w[len] = scal*mb[kb];
                    memset(w[len], 0, nnr*sizeof(FLOAT));
                    memcpy(&w[len][firstrow],&ba[j][k*dimR],dimR*sizeof(FLOAT));//memcpy(csmat->ja[i],ja,nzcount*sizeof(int));
                }
                else{
                    memcpy(&w[jpos][firstrow],&ba[j][k*dimR],dimR*sizeof(FLOAT));//memcpy(csmat->ja[i],ja,nzcount*sizeof(int));//w[jpos] += scal*mb[kb];
                }//if jpos==0


                //mat[bszc[col]+k/dimR][bsz[i]+k%dimR]=ba[j][k];
            }
            //printf("\n w[j]=%f",w[j]); getchar();

        }
        //csmat->ma[i] = (FLOAT*)malloc(nzcount*sizeof(FLOAT));
        //printf("\n csmat->ma[i]=%d",csmat->ma[i]);
        //memcpy(csmat->ma[i],w,nzcount*sizeof(FLOAT));
        //printf("\n after mem"); getchar();
    }

    if (len>=0){
        len++;//important
        mat->ja = (int *)Malloc(sizeof(int)*len,"vb2csr2d mat->ja");//C->ja[ii] = (int *)Malloc(sizeof(int)*len,"csrmulcsr C->ja");
        mat->jab = (int *)Malloc(sizeof(int)*(len+1),"vb2csr2d mat->ja");//C->ja[ii] = (int *)Malloc(sizeof(int)*len,"csrmulcsr C->ja");
        mat->ma =(BData *)Malloc(sizeof(BData)*len,"vb2csr2d mat->ma");
        memcpy(mat->ja, jc, len*sizeof(int));
        memcpy(mat->jab, jcb, len*sizeof(int));
        mat->jab[len] = -1;
        //output_intvector("indexforcolumn.coo" ,jc,0, len);getchar();
        //output_intvector("bindexforcolumn.coo" ,jcb,0, len);getchar();
        memcpy(mat->ma,  w, len*sizeof(BData));
        //output_dblvector("w[0]forcolumn.coo" ,w[0],0, nnr);getchar();
        mat->nzcount = len;
    }
    //printf("\n len=%d",len);
    //free(w);
    //printf("\n n2=%d",csmat->n); getchar();
    free(w);
    free(iw);
    free(jc);
    free(jcb);
    return 0;
}


int csr2d2vb(cs2dptr mat, vbsptr vbmat, int n, int nc, int *bsz, int *bszc)
{
    /*
 * This subroutine converts the sparse 2D FLOAT pointer format to vbsptr matrix
 *
 *----------------------------------------------------------------------
typedef struct Spa2Dmt {
--------------------------------------------- 
| sparse 2D pointer to matrix,
| only store the nonzero columns. 
| specially for F in Barms
|---------------------------------------------
  int n;
  int nc;
  int nzcount;
  int *ja;
  int *jab
  FLOAT **ma;
} Spar2DMat, *cs2dptr;
 * on entry:
 *----------
 *
 * sparse 2D FLOAT pointer format matrix
 *
 * on return:
 *-----------
 *
 *  vbmat = Various Block Sparse Row format Matrix, vbsptr
 *
 *  ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int i, k = 0, nzcount, col,  dimR, dimC, blocksz,ii,kk = 0,*iw,colb;
    int *nB=NULL,*nC=NULL, *jab = mat->jab;
    FLOAT *workbuff=NULL;
    BData *w, *ma = mat->ma;

    nB = (int*)malloc(n*sizeof(int));
    for(i=0;i<n;i++)
    {
        nB[i] = B_DIM(bsz,i);//memset(mat[i], 0, nnr*sizeof(FLOAT));
    }
    nC = (int*)malloc(nc*sizeof(int));
    for(i=0;i<nc;i++)
    {
        nC[i] = B_DIM(bszc,i);//memset(mat[i], 0, nnr*sizeof(FLOAT));
    }
    //printf("\n n=%d,nc=%d",n,nc); getchar();
    //output_intvector("nBin.coo" ,nB,0,n );
    //output_intvector("nCin.coo" ,nC,0,nc );
    if (setupVBMat1(vbmat,n,nc,nB,nC)) return 1;
    //setupCS(csmat,n,1);
    //printf("\n vbmat->bsz=%d",vbmat->bszc[5]); getchar();

    iw =(int*)malloc(nc*sizeof(int));
    w =(BData*)malloc(nc*sizeof(BData));

    for( i = 0; i < n; i++ )//n is the row dimension
    {
        //printf("\n i=%d",i); getchar();
        dimR = B_DIM(bsz,i);//get the row block dimension of the block
        nzcount = vbmat->nzcount[i] = 0;//
        //printf("\n csmat->nzcount[i]=%d",csmat->nzcount[i]); getchar();
        //ja = vbmat->ja[i];
        //printf("\n vbmat->ja[i]=%d",vbmat->ja[i]); getchar();
        //ba = vbmat->ba[i];
        //printf("\n n=%d,nc=%d",n,nc); getchar();
        for( col = 0; col < mat->nzcount; col++ ) //nc is the column block dimension
        {
            //printf("\n j=%d nzcount=%d",j,nzcount); getchar();
            //col = ja[j];
            colb = jab[col];
            dimC = B_DIM(bszc,colb);//obtain the column dimension of the block
            blocksz = dimR*dimC;//caculate the blocksize
            //printf("\n blocksz=%d",blocksz); getchar();
            if(col==0){
                workbuff = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                memset(workbuff, 0, blocksz*sizeof(FLOAT)); //??
                k = 0, kk = 0;
                //printf("\n iw[col]=%d",iw[col]); getchar();
            }
            else if(jab[col]!=jab[col-1]){
                workbuff = (FLOAT*)malloc(blocksz*sizeof(FLOAT));
                memset(workbuff, 0, blocksz*sizeof(FLOAT));
                k = 0, kk = 0;
            }

            //for(k=0,kk=0,jj=0; jj< dimC;jj++)
            //{
            for(ii=0;ii<dimR;ii++){
                //printf("\n i=%d,col=%d,mat[jj][ii]=%f",i,col,mat[jj][ii]); getchar();
                if(ma[col][bsz[i]+ii]){
                    workbuff[k] = ma[col][bsz[i]+ii];kk++;
                }
                k++;
            }
            //mat[bszc[col]+k%dimC][bsz[i]+k/dimC]=ba[j][k];
            //}
            //printf("\n kk=%d,jab[col+1]=%d,jab[col]=%d",kk,jab[col+1],jab[col]); getchar();
            //output_dblvector("workbuff.coo" ,workbuff,0, blocksz);getchar();
            if((jab[col+1]!=jab[col])&&kk){
                iw[nzcount]=colb;w[nzcount]=workbuff;nzcount++;
            }
            else if(jab[col+1]!=jab[col]) free(workbuff);
            //if(kk){//??

            //	iw[nzcount]=colb;w[nzcount]=workbuff;nzcount++;
            //printf("\n iw[col]=%d",iw[col]); getchar();
            //}
            //else free(workbuff);
            //}

        }
        //printf("\n nzcount=%d",nzcount); getchar();
        vbmat->nzcount[i]=nzcount;
        if( nzcount > 0 ){
            vbmat->ja[i] =(int*)malloc(nzcount*sizeof(int));
            vbmat->ba[i] =(BData*)malloc(nzcount*sizeof(BData));
            memcpy(vbmat->ja[i],iw,nzcount*sizeof(int));
            //printf("\n vbmat->ja[i]=%d",vbmat->ja[i][0]);
            memcpy(vbmat->ba[i],w,nzcount*sizeof(BData));//memcpy(Aflit->ja[i],iw,len*sizeof(int));
            //csmat->ma[i] = (FLOAT*)malloc(nzcount*sizeof(FLOAT));
            //printf("\n csmat->ma[i]=%d",csmat->ma[i]);
            //memcpy(csmat->ma[i],w,nzcount*sizeof(FLOAT));
            //printf("\n after mem"); getchar();
        }
    }

    free(w);
    free(iw);
    free(nB);
    free(nC);
    //printf("\n n2=%d",csmat->n); getchar();
    return 0;
}


int cleanCS2D(cs2dptr amat)
{
    /*----------------------------------------------------------------------
| Free up memory allocated for Spa2Dmt structs.
|----------------------------------------------------------------------
typedef struct Spa2Dmt {
--------------------------------------------- 
| sparse 2D pointer to matrix,
| only store the nonzero columns. 
| specially for F in Barms
|---------------------------------------------
  int n;	// row dimension of matrix
  int nc;	 column dimension of matrix
  int nzcount;   length of each row
  int *ja;      pointer to store column indices, sorted
  int *jab;	 pointer to store block column indices, sorted
  FLOAT **ma;   pointer to store full dense column pointer
} Spar2DMat, *cs2dptr;
| on entry:
|==========
| ( amat )  =  Pointer to a SpaFmt struct.
|--------------------------------------------------------------------*/

    int i;
    if (amat == NULL) return 0;
    if (amat->n < 1) return 0;
    //printf("\n cleanCS n=%d",amat->n);// getchar();
    //for (i=0; i<amat->nzcount; i++) {
    //if((i>=(amat->n-3))||(i % 500==0)) { printf("\n cleanCS 1:  row %d",i); getchar(); }
    //printf("\n cleanCS 2");
    //  if( amat->ma[j] ) free(amat->ma[i]);amat->ma[i]=NULL;
    //printf("\n cleanCS before ja[%d]",i);
    //if (amat->ja[i]==NULL) {printf("==NULL,cann't free");  getchar();}
    //free(amat->ja[i]);
    //printf("\n cleanCS after ja[%d]",i); getchar();
    //}
    if( amat->ma ) {
        for( i=0; i<amat->nzcount; i++ ) {
            if( amat->ma[i] ) free( amat->ma[i] );
        }
        free( amat->ma );amat->ma = NULL;
    }

    //printf("\n cleanCS 4");getchar();
    //if (amat->ma) free(amat->ma);amat->ma = NULL;
    //printf("\n cleanCS 5");getchar();
    free(amat->ja);amat->ja = NULL;
    //printf("\n cleanCS 6");getchar();
    free(amat->jab);amat->jab = NULL;
    //printf("\n cleanCS 7");getchar();
    free(amat);
    //printf("\n cleanCS 8"); getchar();
    amat = NULL;//add
    return 0;
}
/*---------------------------------------------------------------------
|     end of cleanCS
|--------------------------------------------------------------------*/



int computLUFsparse(vbsptr F, vbiluptr lu, vbsptr FF){
    //comput the matrix product (LU)^-1*F

    int i, n = F->n, nc = F->nc;
    int *bsz = F->bsz,*bszc = F->bszc,nnr = bsz[n],nzcount;//nnc = bszc[nc],
    FLOAT *workbuff;
    cs2dptr Bmat;//BData *Bmat;
    //~ double elapsedtime;
    //~ time_t tm1,tm2;
    //~ double tm11=0.0, tm22=0.0, tmpp=0.0;
    BData *ma;// = mat->ma;

    Bmat = (cs2dptr)Malloc(sizeof(Spar2DMat),"computLUF");//Bmat = (BData *)Malloc(nnc*sizeof(BData),"computLUF");

    //~ tm1 = time (NULL);
    //~ tm11 = sys_timer();
    vb2csr2d(F, Bmat);//int vb2csr2d(vbsptr vbmat, cs2dptr mat)//ierr = vb2BData(F, Bmat);//int vb2BData(vbsptr vbmat, BData *mat)
    //outputcs2dptr(Bmat,"Bmat.coo",1);
    //~ tm2 = time (NULL);
    //~ tm22 = sys_timer();
    //~ tmpp = tm22 - tm11;
    //~ elapsedtime = difftime (tm2, tm1);
    //printf("\n         cputimevb->||||=%f,calendartimevb->||||=%f\n",tmpp,elapsedtime);


    ma = Bmat->ma;nzcount = Bmat->nzcount;
    //printf("\n         nzcount=%d\n",nzcount);
    workbuff = (FLOAT *)Malloc(nnr*sizeof(FLOAT),"computLUF");

    //~ tm1 = time (NULL);
    //~ tm11 = sys_timer();

    //outputcs2dptr(Bmat,"Bmat[i].coo",1);

    for( i = 0; i < nzcount; i++ )
    {
        //output_dblvector("ma[i].coo" ,ma[i],0, nnr);getchar();
        vblusolC(ma[i],workbuff,lu);//int vblusolC( FLOAT *y, FLOAT *x, vbiluptr lu)
        copyBData(nnr,1,ma[i],workbuff,0);//void copyBData( int m, int n, BData dst, BData src, int isig )
        //output_dblvector("ma[i]aftersolve.coo" ,ma[i],0, nnr);getchar();
    }
    //outputcs2dptr(Bmat,"Bmataftersolve.coo",1);
    //~ tm2 = time (NULL);
    //~ tm22 = sys_timer();
    //~ tmpp = tm22 - tm11;
    //~ elapsedtime = difftime (tm2, tm1);
    // printf("\n         cputimeB^(-1)*|=%f,calendartimeB^(-1)*|=%f\n",tmpp,elapsedtime);


    /*nB = (int*)malloc(n*sizeof(int));
  for(i=0;i<n;i++){
      nB[i] = B_DIM(bsz,i);//memset(mat[i], 0, nnr*sizeof(FLOAT));
  }
  nC = (int*)malloc(nc*sizeof(int));
  for(i=0;i<nc;i++){
      nC[i] = B_DIM(bszc,i);//memset(mat[i], 0, nnr*sizeof(FLOAT));
  }
  if (setupVBMat1(FF,n,nc,nB,nC)) return 1;*/
    //outputBData( Bmat, "BFmat.coo", 1,nnr,nnc );//int outputBData( BData *mat, char *filename, int onebase,int n,int nc)
    //~ tm1 = time (NULL);
    //~ tm11 = sys_timer();
    //printf("\n         cputimeB^(-1)*|=%f,calendartimeB^(-1)*|=%f\n",tmpp,elapsedtime);
    //output_intvector("ja.coo" ,Bmat->ja,0, Bmat->nzcount);getchar();
    //output_intvector("jab.coo" ,Bmat->jab,0, Bmat->nzcount);getchar();
    qqsortm(Bmat->ja,Bmat->jab, Bmat->ma, 0, Bmat->nzcount-1);//qqsortm(Fs2D->ja,Fs2D->jab, Fs2D->ma, 0, Fs2D->nzcount);
    //output_intvector("jasorted.coo" ,Bmat->ja,0, Bmat->nzcount);getchar();
    //output_intvector("jabsorted.coo" ,Bmat->jab,0, Bmat->nzcount);getchar();
    csr2d2vb(Bmat, FF, n, nc, bsz, bszc);//int csr2d2vb(cs2dptr mat, vbsptr vbmat, int n, int nc, int *bsz, int *bszc)//ierr = BData2vb(Bmat, FF, n, nc, bsz, bszc);//int BData2vb(BData *mat, vbsptr vbmat, int n, int nc, int *bsz, int *bszc)
    //~ tm2 = time (NULL);
    //~ tm22 = sys_timer();
    //~ tmpp = tm22 - tm11;
    //~ elapsedtime = difftime (tm2, tm1);
    //printf("\n         cputime||||->vb=%f,calendartime||||->vb=%f\n",tmpp,elapsedtime);


    //for( i = 0; i < nnc; i++ )
    //{
    //  free(Bmat[i]);
    //}
    //free(Bmat);
    cleanCS2D(Bmat);//int cleanCS2D(cs2dptr amat)
    free(workbuff);
    //printf("\n n2=%d",csmat->n); getchar();
    //free(nB);
    //free(nC);
    return 0;
}

int computLUFinpartition( vbbsptr Bb, vbbsptr Ff, vbbsptr FF, vbbiluptr luu,double *droptol, int *lfil, FILE *flog)
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
    int i,j,max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(FLOAT);
    int nset,*nBB;
    double tol;
    BData *w;
    //vbbiluptr luu = NULL;
    vbsptr *vba,*vbaF,*vbaFF;
    //~ double elapsedtime;
    //~ time_t tm1,tm2;
    //~ double tm11=0.0, tm22=0.0, tmpp=0.0;
    //vbsptr FF = NULL,FFF = NULL;
    //double elapsedtime;
    //time_t tm1,tm2;
    //double tm11=0.0, tm22=0.0, tmpp=0.0;

    /*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/
    //tm1 = time (NULL);
    //tm11 = sys_timer();
    //printf("\n  computLUFinpartition 11");getchar();
    tol  = droptol[5];
    nset = Bb->n;
    //printf("\n  nset=%d",nset);getchar();
    //luu = (vbbiluptr)Malloc( sizeof(VBBILUSpar), "main" );

    nBB = (int *) Malloc(nset*sizeof(int), "vbSplit4copy:0" );
    //output_intvector("bszinFF.coo" ,Ff->bsz,0, nset+1);getchar();
    for(i=0;i<nset;i++)
    {
        nBB[i] = B_DIM(Ff->bsz,i);//memset(mat[i], 0, nnr*sizeof(double));
    }
    //printf("\n  computLUFinpartition 22");getchar();
    //output_intvector("nBBinFF.coo" ,nBB,0, nset);getchar();

    if (setupVBvector(FF,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
    if (setupVBilutvector(luu,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
    //if (setupVBvector(FF,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
    //FF->n = nset;luu->n = nset;
    //~ tm1 = time (NULL);
    //~ tm11 = sys_timer();
    for( i = 0; i < nset; i++ ){
        //printf("\n  i=%d",i);getchar();
        vba = Bb->vba;
        vbaF = Ff->vba;
        w = (BData *)Malloc(vba[i]->n*sizeof(BData),"main");
        for( j = 0; j < vba[i]->n; j++ )
            w[j] = (FLOAT *)Malloc( max_blk_sz, "main" );
        //printf("\n  luu->vlu=%d",luu->vlu);getchar();
        luu->vlu[i] = (vbiluptr)Malloc( sizeof(VBILUSpar), "main" );
        //printf("\n  i=%d",i);getchar();
        //outputvbmat(vba[i],"vba[i].coo",1);
        vbilutC( vba[i], luu->vlu[i], lfil[5], tol, w, flog );//int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp )
        //ierr = vbilukC( fill_lev, vba[i], luu->vlu[i], flog );//vbchange
        //	    FF = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
        vbaFF = FF->vba;
        //outputvbmat1(vbaF[i],"vbaF[i].coo",1);
        //printf("\n  vbaFF[i]=%d",vbaFF[i]);getchar();
        vbaFF[i] = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
        computLUFsparse(vbaF[i],luu->vlu[i],vbaFF[i]);//ierr = computLUFsparse(vbaF[i],luu->vlu[i],vbaFF[i]);//int computLUF(vbsptr F, vbiluptr lu, vbsptr FF)
        //printf("\n  vbaFF[i]=%d",vbaFF[i]);getchar();
        //outputvbmat1(vbaFF[i],"vbaFF[i].coo",1);
        for( j = 0; j < vba[i]->n; j++ )
            free( w[j] );
        free( w );
    }
    //~ tm2 = time (NULL);
    //~ tm22 = sys_timer();
    //~ tmpp = tm22 - tm11;
    //~ elapsedtime = difftime (tm2, tm1);
    //printf("\n     cputimeforincomputLUF=%f,calendartimeforincomputLUF=%f\n",tmpp,elapsedtime);
    
    if (nBB) free(nBB);
    //cleanVBMat1( FF );
    //cleanVBMat( FFF );
    return 0;
}

int Ffassemble( vbsptr FF, vbbsptr Ff)
{

    /*typedef struct VBBSpaFmt {
    int n;
    int *bsz;
    vbsptr *vba;
} VBBSparMat, *vbbsptr;*/
    /*typedef struct VBSpaFmt {
    int n;
    int *bsz;
    int *nzcount;
    int **ja;
    BData **ba;
    BData *D;
} VBSparMat, *vbsptr;*/     

    int j,i,k,*ja, *nB, *nC, dimR, dimC,firstFi,rowFF,nzcount, col,*FFja,nset = Ff->n;
    int *bszz = Ff->bsz,n = bszz[nset],nc = Ff->vba[0]->nc,*bszc = Ff->vba[0]->bszc,*bsz;
    BData *ba,*FFba;
    vbsptr *vba;
    vbsptr Fi;
    // printf("\nFf->n=%d,n=%d",Ff->n,n); getchar();

    nB = (int *) Malloc(n*sizeof(int), "vbSplit4copy:0" );
    nC = (int *) Malloc(nc*sizeof(int), "vbSplit4copy:0" );

    for (i=0; i<nset; i++){
        //nB is size of diagonal fine block, B->n is the dimension of it.
        vba = Ff->vba;
        for (j=0; j<B_DIM(bszz,i); j++){
            //printf("\ni=%d,j=%d,B_DIM(bszz,i)=%d",i,j,B_DIM(bszz,i)); getchar();
            nB[bszz[i]+j]=B_DIM(vba[i]->bsz,j);
        }
    }
    for (i=0; i<nc; i++)
        nC[i]=B_DIM(bszc,i);
    //output_intvector("nBinFfassemble.coo" ,nB,0, n+1);getchar();
    if (setupVBMat1(FF,n,nc,nB,nC)) return 1; //int setupVBMat1( vbsptr vbmat, int n,int nc, int *nBR, int *nBC )
    //nBB is size of diagonal coarse block, nset is the dimension of it.
    //if (setupVBvector(Ff,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
    //free(nB);
    //free(nC);
    //printf("\nFF->nc=%d,FF->n=%d",FF->nc,FF->n); getchar();
    //output_intvector("FF->bsz.coo" ,FF->bsz,0, FF->n+1);getchar();
    // output_intvector("FF->bszc.coo" ,FF->bszc,0, FF->nc+1);
    //printf("\nafter setup "); getchar();
    bsz = FF->bsz;
    for( i = 0; i < nset; i++ ) {
        //printf("\ni=%d",i); getchar();
        vba = Ff->vba;
        Fi = vba[i];
        firstFi = bszz[i];//in the fine block vector, find the first address of each coarse block

        //printf("\nfirstBi=%d",firstFi); getchar();
        //if (setupVBMat1(F,bsize,csize,nB,nC)) goto label111;
        //output_intvector("bszinBi.coo" ,Bi->bsz,0, nBB[i]+1);getchar();
        //printf("\Bi->n=%d",Bi->n); getchar();

        for( j = 0; j < B_DIM(bszz,i); j++ ) {
            rowFF = firstFi+j;
            dimR = B_DIM(bsz,rowFF);
            nzcount = FF->nzcount[rowFF] = Fi->nzcount[j];
            ja = Fi->ja[j];
            ba = Fi->ba[j];
            if (nzcount>0) {
                FF->ja[rowFF] = (int *) Malloc(nzcount*sizeof(int), "vbSplit4copy:5" );
                FF->ba[rowFF] = (BData *) Malloc(nzcount*sizeof(BData), "vbSplit4copy:6" );
            }
            //printf("\ni=%d,j=%d,nz=%d",i,j,nzcount); getchar();
            // memcpy(Bi->ja[j], ja, nzcount*sizeof(int));//ja-
            //printf("\ni=%d,j=%d",i,j); getchar();

            FFja = FF->ja[rowFF];
            FFba = FF->ba[rowFF];
            for( k = 0; k < nzcount; k++ ) {

                col = ja[k];
                //printf("\ni=%d,j=%d,rowFF=%d,k=%d,col=%d",i,j,rowFF,k,col); getchar();
                dimC = B_DIM(bszc,col);
                FFja[k] = ja[k];
                FFba[k] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
                copyBData( dimC, dimR, FFba[k], ba[k], 0 );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
            }
        }
        //printf("\n vba[i]->n=%d",vba[i]->n); getchar();
    }
    // Bb->maxsize = maxsize;
    //Ff->maxsize = maxsize;
    // printf("\n Bb->maxsize=%d",Bb->maxsize); getchar();
    if (nB) free(nB);
    if (nC) free(nC);
    return 0;
}


int cleanVBBMat( vbbsptr vbbmat ){
    /*----------------------------------------------------------------------
| Free up memory allocated for VBILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|   ( vbbmat )  =  Pointer to a VBILUSpar struct.
|--------------------------------------------------------------------*/

    /*typedef struct VBBSpaFmt {
    int n;
    int maxsize;
    int *bsz;
    vbsptr *vba;
} VBBSparMat, *vbbsptr;*/

    /*  if (amat->lu) {
    cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    amat->lu = NULL;*/

    int n = vbbmat->n, i;
    if( NULL == vbbmat ) return 0;
    if( vbbmat->vba ) {
        for( i = 0; i < n; i++ ) {
            if( vbbmat->vba[i] ) {
                cleanVBMat1(vbbmat->vba[i]);
                vbbmat->vba[i]=NULL;
            }
        }
        free( vbbmat->vba );
    }
    if( vbbmat->bsz ) free( vbbmat->bsz );
    free( vbbmat );
    return 0;
}


int cleanVBBILU( vbbiluptr luu ){
    /*----------------------------------------------------------------------
| Free up memory allocated for VBILUSpar structs.
|----------------------------------------------------------------------
| on entry:
|   ( luu )  =  Pointer to a VBILUSpar struct.
|--------------------------------------------------------------------*/

    /*typedef struct VBBILUfac {
    int n;
    int *bsz;
    vbiluptr *vlu;
} VBBILUSpar, *vbbiluptr;*/

    /*  if (amat->lu) {
    cleanVBILU(amat->lu);//cleanVBILU(vbiluptr lu);
    amat->lu = NULL;*/

    int n = luu->n, i;
    if( NULL == luu ) return 0;
    if( luu->vlu ) {
        for( i = 0; i < n; i++ ) {
            if( luu->vlu[i] ) { cleanVBILU(luu->vlu[i]);luu->vlu[i]=NULL; }
        }
        free( luu->vlu );
    }
    if( luu->bsz ) free( luu->bsz );
    free( luu );
    return 0;
}

int BLUassemble( vbiluptr lu, vbbiluptr luu)//int Ffassemble( vbsptr FF, vbbsptr Ff)
{

    /*typedef struct VBBILUfac {
    int n;
    int *bsz;
    vbiluptr *vlu;
} VBBILUSpar, *vbbiluptr;*/

    int j,i,k,*Lija,*Lja,*Uija,*Uja, *nB, dimR, dimC,firstFi,rowlu,nzcountL,nzcountU, col,nset = luu->n,*bszforlu;
    int *bszz = luu->bsz,n = bszz[nset],*bsz;
    BData *Liba,*Uiba,*Lba,*Uba,*D,*Di;
    vbiluptr *vlu;
    vbiluptr lui;
    vbsptr L,U,Li,Ui;
    //printf("\nFf->n=%d,n=%d",luu->n,n); getchar();

    nB = (int *) Malloc(n*sizeof(int), "vbSplit4copy:0" );
    //nC = (int *) Malloc(nc*sizeof(int), "vbSplit4copy:0" );

    for (i=0; i<nset; i++){
        //nB is size of diagonal fine block, B->n is the dimension of it.
        vlu = luu->vlu;
        for (j=0; j<B_DIM(bszz,i); j++){
            //printf("\ni=%d,j=%d,B_DIM(bszz,i)=%d",i,j,B_DIM(bszz,i)); getchar();
            nB[bszz[i]+j]=B_DIM(vlu[i]->bsz,j);
        }
    }

    if( nB ) {
        bszforlu = (int *)Malloc( sizeof(int)*(n+1), "setupVBMat" );
        bszforlu[0] = 0;
        for( i = 1; i <= n; i++ ) {
            bszforlu[i] = bszforlu[i-1] + nB[i-1];
        }
    } else
        bszforlu = NULL;

    //output_intvector("nBinFfassemble.coo" ,nB,0, n+1);getchar();
    setupVBILU( lu, n, bszforlu );
    L = lu->L;
    U = lu->U;
    D = lu->D;  //if (setupVBMat1(FF,n,nc,nB,nC)) return 1; //int setupVBMat1( vbsptr vbmat, int n,int nc, int *nBR, int *nBC )
    //nBB is size of diagonal coarse block, nset is the dimension of it.
    //printf("\nFF->nc=%d,lu->n=%d",lu->nc,lu->n); getchar();
    //output_intvector("lu->bsz.coo" ,lu->bsz,0, lu->n+1);getchar();
    //printf("\nafter setup "); getchar();
    bsz = lu->bsz;
    for( i = 0; i < nset; i++ ) {
        //printf("\ni=%d",i); getchar();
        vlu = luu->vlu;
        lui = vlu[i];
        Li = lui->L;
        Ui = lui->U;
        Di = lui->D;
        firstFi = bszz[i];//in the fine block vector, find the first address of each coarse block

        //printf("\nfirstBi=%d",firstFi); getchar();

        //if (setupVBMat1(F,bsize,csize,nB,nC)) goto label111;
        //output_intvector("bszinBi.coo" ,Bi->bsz,0, nBB[i]+1);getchar();
        //printf("\Bi->n=%d",Bi->n); getchar();

        for( j = 0; j < B_DIM(bszz,i); j++ ) {
            rowlu = firstFi+j;
            dimR = B_DIM(bsz,rowlu);
            D[rowlu] = (BData) Malloc(dimR*dimR*sizeof(FLOAT), "vbSplit4copy:6" );
            //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
            copyBData( dimR, dimR, D[rowlu], Di[j], 0 );
            nzcountL = L->nzcount[rowlu] = Li->nzcount[j];
            nzcountU = U->nzcount[rowlu] = Ui->nzcount[j];
            Lija = Li->ja[j];
            Liba = Li->ba[j];
            Uija = Ui->ja[j];
            Uiba = Ui->ba[j];
            if (nzcountL>0) {
                L->ja[rowlu] = (int *) Malloc(nzcountL*sizeof(int), "vbSplit4copy:5" );
                L->ba[rowlu] = (BData *) Malloc(nzcountL*sizeof(BData), "vbSplit4copy:6" );
            }
            if (nzcountU>0) {
                U->ja[rowlu] = (int *) Malloc(nzcountU*sizeof(int), "vbSplit4copy:5" );
                U->ba[rowlu] = (BData *) Malloc(nzcountU*sizeof(BData), "vbSplit4copy:6" );
            }
            //printf("\ni=%d,j=%d,nz=%d",i,j,nzcount); getchar();
            // memcpy(Bi->ja[j], ja, nzcount*sizeof(int));//ja-
            //printf("\ni=%d,j=%d",i,j); getchar();

            Lja = L->ja[rowlu];
            Lba = L->ba[rowlu];
            Uja = U->ja[rowlu];
            Uba = U->ba[rowlu];
            for( k = 0; k < nzcountL; k++ ) {

                col = Lija[k]+firstFi;
                //printf("\ni=%d,j=%d,rowlu=%d,k=%d,col=%d",i,j,rowlu,k,col); getchar();
                dimC = B_DIM(bsz,col);
                Lja[k] = Lija[k]+firstFi;
                Lba[k] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
                copyBData( dimC, dimR, Lba[k], Liba[k], 0 );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
            }
            for( k = 0; k < nzcountU; k++ ) {

                col = Uija[k]+firstFi;
                //printf("\ni=%d,j=%d,rowlu=%d,k=%d,col=%d",i,j,rowlu,k,col); getchar();
                dimC = B_DIM(bsz,col);
                Uja[k] = Uija[k]+firstFi;
                Uba[k] = (BData) Malloc(dimR*dimC*sizeof(FLOAT), "vbSplit4copy:6" );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
                copyBData( dimC, dimR, Uba[k], Uiba[k], 0 );
                //printf("\ni=%d,j=%d,k=%d",i,j,k); getchar();
            }
        }
        //printf("\n vlu[i]->n=%d",vlu[i]->n); getchar();
    }
    // Bb->maxsize = maxsize;
    //Ff->maxsize = maxsize;
    // printf("\n Bb->maxsize=%d",Bb->maxsize); getchar();
    lu->DiagOpt = 1;//vbchange
    if (nB) free(nB);
    if (bszforlu) free(bszforlu);
    return 0;
}


int computschurpartition(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol, int *lfil, vbsptr schur, int *nBB, int nset)
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
    //,max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(FLOAT);
    //int lfil;
    double tol = droptol[5];
    // BData *w;
    vbiluptr lu = NULL;
    vbsptr FF = NULL,FFF = NULL;
    vbbsptr Ff,Fff,Bb;
    vbbiluptr luu;
    //~ double elapsedtime;
    //~ time_t tm1,tm2;
    //~ double tm11=0.0, tm22=0.0, tmpp=0.0;
    /*---------------------------------------------------------------------
|     Sort the matrix and separate into   |  B  F  |
|                                         |        |
|                                         |  E  C  |
|--------------------------------------------------------------------*/

    //w = (BData *)Malloc(B->n*sizeof(BData),"main");
    //for( i = 0; i < B->n; i++ )
    //w[i] = (FLOAT *)Malloc( max_blk_sz, "main" );
    //lfil = io.lfil0;
    //printf("\n lfil=%d,tol=%20.16e",lfil,tol); getchar();

    // ierr = vbilutC( B, lu, lfil, tol, w, flog );//int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp )
    //ierr = setupVBP4 (vbmat, B->n, C->n, F, E , int *bsz,vbiluptr lu );//int setupVBP4 (vbp4ptr vbmat, int Bn, int Cn,  vbsptr F,  vbsptr E , int *bsz,vbiluptr lu )
    FF = (vbsptr) Malloc(sizeof(VBSparMat), "arms2:7" );
    Ff = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
    Fff = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
    Bb = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
    luu = (vbbiluptr) Malloc(sizeof(VBBILUSpar), "arms2:7" );
    lu = (vbiluptr)Malloc( sizeof(VBILUSpar), "main" );
    /*extern int Bsplit(vbsptr B, vbbsptr Bb, vbsptr F, vbbsptr Ff, int *nBB, int nset);
extern int computLUFinpartition( vbbsptr Bb, vbbsptr Ff, vbbsptr FF, vbbiluptr luu,io_t io, FILE *flog);
extern int Ffassemble( vbsptr FF, vbbsptr Ff);
extern int BLUassemble( vbiluptr lu, vbbiluptr luu);*/
    //~ tm1 = time (NULL);
    //~ tm11 = sys_timer();
    Bsplit(B, Bb, vbmat->F, Ff, nBB, nset);
    //~ tm2 = time (NULL);
    //~ tm22 = sys_timer();
    //~ tmpp = tm22 - tm11;
    //~ elapsedtime = difftime (tm2, tm1);
    //printf("\n     cputimeBsplit=%f,calendartimeBsplit=%f\n",tmpp,elapsedtime);
    //~ tm1 = time (NULL);
    //~ tm11 = sys_timer();
    computLUFinpartition( Bb, Ff, Fff, luu, droptol, lfil, stderr);
    //~ tm2 = time (NULL);
    //~ tm22 = sys_timer();
    //~ tmpp = tm22 - tm11;
    //~ elapsedtime = difftime (tm2, tm1);
    //printf("\n     cputimecomputLUF=%f,calendartimecomputLUF=%f\n",tmpp,elapsedtime);
    cleanVBBMat(Bb);
    cleanVBBMat(Ff);
    Ffassemble( FF, Fff);
    cleanVBBMat(Fff);//int cleanVBBMat( vbbsptr vbbmat )
    BLUassemble( lu, luu);
    cleanVBBILU(luu);//int cleanVBBILU( vbbiluptr luu )
    vbmat->lu = lu;

    //ierr = computLUF(vbmat->F,vbmat->lu,FF);//int computLUF(vbsptr F, vbiluptr lu, vbsptr FF)

    //printf("\n     cputimeB^(-1)*F=%f,calendartimeB^(-1)*F=%f\n",tmpp,elapsedtime);
    //outputvbmat1(FF,"FFinB.coo",1);
    FFF = (vbsptr) Malloc(sizeof(VBSparMat), "arms2:7" );
    //outputvbmat1(vbmat->E,"vbmat->E.coo",1);

    vbmulvb(1,vbmat->E,FF,FFF);//int vbmulvb(int job, vbsptr A, vbsptr B, vbsptr C)
    vbplussvb1(1,C,-1.0,FFF,schur,tol);//int vbplussvb(int job, vbsptr A, FLOAT s, vbsptr B, vbsptr C)??
    
    //outputvbmat(schur,"schur.coo",1);
    
    //for( i = 0; i < B->n; i++ )
    //free( w[i] );
    //free( w );
    cleanVBMat1( FF );
    cleanVBMat( FFF );
    return 0;
}
//~ 
//~ int computLUFinpartitionp( vbbsptr Bb, vbbsptr Ff, vbbsptr FF, iilutptr luu, double *droptol, int *lfil, int *ipar, FILE *flog)
//~ {
//~ /*---------------------------------------------------------------------
//~ | Convert permuted vbspaFmt struct to VBPerMat4 struct 
//~ |                - matrix already permuted
//~ |----------------------------------------------------------------------
//~ | on entry:
//~ |========== 
//~ | ( vbmat )  =  Matrix stored in vbspaFmt format.
//~ |              Internal pointers (and associated memory) destroyed before
//~ |              return.
//~ |
//~ | On return:
//~ |===========
//~ |
//~ | B, E, F, C = 4 blocks in 
//~ | 
//~ |          | B   F |      
//~ |   vbmat= |       | 
//~ |          | E   C | 
//~ | 
//~ |
//~ |       integer value returned:
//~ |             0   --> successful return.
//~ |             1   --> memory allocation error.
//~ |--------------------------------------------------------------------*/
//~ int i,j,ierr = 0;//max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(FLOAT);
//~ int nset,*nBB,fill_lev,nbic,pivoting = 1;
//~ double tol,PERMTOL = 0.99;
//~ //BData *w;
//~ csptr Bic;
//~ //vbbiluptr luu = NULL;
//~ vbsptr *vba,*vbaF,*vbaFF;
//~ ilutptr ilsch;
//~ //vbsptr FF = NULL,FFF = NULL;
//~ //double elapsedtime;
//~ //time_t tm1,tm2;
//~ //double tm11=0.0, tm22=0.0, tmpp=0.0;
//~
//~ /*---------------------------------------------------------------------
//~ |     Sort the matrix and separate into   |  B  F  |
//~ |                                         |        |
//~ |                                         |  E  C  |
//~ |--------------------------------------------------------------------*/
//~ //tm1 = time (NULL);
//~ //tm11 = sys_timer();
//~ //printf("\n  computLUFinpartition 11");getchar();
//~ //printf("\n  diagscal=%d",diagscal);getchar();
//~ //io.diagscal = diagscal;
//~ //lfil = io.lfil0;
//~ //tol  = io.tol0;
//~ nset = Bb->n;
//~ //printf("\n  nset=%d",nset);getchar();
//~ //luu = (vbbiluptr)Malloc( sizeof(VBBILUSpar), "main" );
//~ 
//~ nBB = (int *) Malloc(nset*sizeof(int), "vbSplit4copy:0" );
//~ //output_intvector("bszinFF.coo" ,Ff->bsz,0, nset+1);getchar();
//~ for(i=0;i<nset;i++)
//~ {
//~ nBB[i] = B_DIM(Ff->bsz,i);//memset(mat[i], 0, nnr*sizeof(FLOAT));
//~ }
//~ //printf("\n  computLUFinpartition 22");getchar();
//~ //output_intvector("nBBinFF.coo" ,nBB,0, nset);getchar();
//~ 
//~ if (setupVBvector(FF,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
//~ if (setupilutpvector(luu,nset,nBB)) return 1;//if (setupVBilutvector(luu,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
//~ //if (setupVBvector(FF,nset,nBB)) return 1;//int setupVBvector( vbbsptr Bb, int nset, int *nBB)
//~ //FF->n = nset;luu->n = nset;
//~ for( i = 0; i < nset; i++ ){
//~ //printf("\n  i=%d",i);getchar();
//~ vba = Bb->vba;
//~ vbaF = Ff->vba;
//~ Bic = (csptr)Malloc( sizeof(SparMat), "main" );
//~ ierr = vbsrc2csr1(vba[i], Bic);
//~ //outputcsmat(Bic,"Bic.coo",1);
//~ nbic = Bic->n;
//~ luu->llu[i] = (ilutptr)Malloc( sizeof(IluSpar), "main" );
//~ setupILUT(luu->llu[i],nbic);
//~ ilsch = luu->llu[i];
//~ //w = (BData *)Malloc(vba[i]->n*sizeof(BData),"main");
//~ //for( j = 0; j < vba[i]->n; j++ )
//~ //w[j] = (FLOAT *)Malloc( max_blk_sz, "main" );
//~ //printf("\n  luu->vlu=%d",luu->vlu);getchar();
//~ ilsch->D1 = NULL;
//~ ilsch->D2 = NULL;
//~ ilsch->perm2 = NULL;
//~ ilsch->rperm = NULL;
//~ ilsch->perm  = NULL;
//~ 
//~ 
//~ if (unsyorder) {
//~ //printf("\n  nbic=%d",nbic);getchar();
//~ int nB,bsize = io.Bsizeinset;
//~ //printf("\n  bsize=%d\n",bsize);getchar();
//~ ilsch->rperm = (int *) Malloc(nbic*sizeof(int), "arms2:ilutpC" );
//~ ilsch->perm = (int *) Malloc(nbic*sizeof(int), "arms2:ilutpC" );
//~ double tolind = 0.0;
//~ PQperm(Bic, 10, ilsch->rperm, ilsch->perm, &nB, tolind) ;
//~ rpermC(Bic,ilsch->rperm);
//~ cpermC(Bic,ilsch->perm);
//~ }
//~ 
//~ 
//~ if (pivoting == 0)
//~ ierr = ilutD(Bic, droptol, lfil, ilsch);
//~ else {
//~ ilsch->perm2 = (int *) Malloc(nbic*sizeof(int), "arms2:ilutpC" );
//~ for (j=0; j<nbic; j++)
//~ ilsch->perm2[j] = j;
//~ //printf("\n  nbic=%d",nbic);getchar();
//~ ierr = ilutpC(Bic, droptol, lfil, PERMTOL,0 , ilsch);//ilutpchange//nbic
//~ }
//~ cleanCS(Bic);//int cleanCS(csptr amat)
//~ //outputcsmat(ilsch->L,"ilsch->L.coo",1);
//~ //outputcsmat(ilsch->U,"ilsch->U.coo",1);
//~ //output_intvector("perm2.coo" ,ilsch->perm2,0, nbic);getchar();
//~ //ierr = arms2(csmat, ipar, droptol, lfil_arr, tolind, ArmsSt, flog);
//~ //printf("\n  i=%d",i);getchar();	 
//~ //outputvbmat(vba[i],"vba[i].coo",1);
//~
//~ /*   setupILUT(ilsch,nC); 
//~ if (ilev > 0) ilsch->C=C;
//~ ilsch->perm2 = NULL;
//~ if (methS[1] == 0)
//~ ierr = ilutD(schur, droptol, lfil, ilsch);
//~ else {
//~ ilsch->perm2 = (int *) Malloc(nC*sizeof(int), "arms2:ilutpC" );
//~ for (j=0; j<nC; j++)
//~ ilsch->perm2[j] = j;
//~ ierr = ilutpC(schur, droptol, lfil, PERMTOL, nC, ilsch);
//~ }*/
//~ //int ilutpC(csptr amat, double *droptol, int *lfil, double permtol,int mband, ilutptr ilusch)
//~ //ierr = vbilutC( vba[i], luu->vlu[i], lfil, tol, w, flog );//int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp )
//~ //ierr = vbilukC( fill_lev, vba[i], luu->vlu[i], flog );//vbchange
//~ //	    FF = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
//~ vbaFF = FF->vba;
//~ //outputvbmat1(vbaF[i],"vbaF[i].coo",1);
//~ //printf("\n  vbaFF[i]=%d",vbaFF[i]);getchar();
//~ vbaFF[i] = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
//~ //printf("\n  vbaFF[i]=%d",vbaFF[i]);getchar();
//~ ierr = computLUFpsparse(vbaF[i],ilsch,vbaFF[i]);//int computLUF(vbsptr F, vbiluptr lu, vbsptr FF)
//~ //printf("\n  vbaFF[i]=%d",vbaFF[i]);getchar();
//~ //outputvbmat1(vbaFF[i],"vbaFF[i].coo",1);
//~ //for( j = 0; j < vba[i]->n; j++ )
//~ //free( w[j] );
//~ //free( w );
//~ }
//~
//~ if (nBB) free(nBB);
//~ //cleanVBMat1( FF );
//~ //cleanVBMat( FFF );
//~ return 0;
//~ }
//~ 
//~ 
//~ int computschurpartitionp(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol, int *lfil, int *ipar, vbsptr schur, int *nBB, int nset)
//~ {
//~ /*---------------------------------------------------------------------
//~ | Convert permuted vbspaFmt struct to VBPerMat4 struct 
//~ |                - matrix already permuted
//~ |----------------------------------------------------------------------
//~ | on entry:
//~ |========== 
//~ | ( vbmat )  =  Matrix stored in vbspaFmt format.
//~ |              Internal pointers (and associated memory) destroyed before
//~ |              return.
//~ |
//~ | On return:
//~ |===========
//~ |
//~ | B, E, F, C = 4 blocks in 
//~ | 
//~ |          | B   F |      
//~ |   vbmat= |       | 
//~ |          | E   C | 
//~ | 
//~ |
//~ |       integer value returned:
//~ |             0   --> successful return.
//~ |             1   --> memory allocation error.
//~ |--------------------------------------------------------------------*/
//~ int ierr = 0;//,max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(FLOAT);
//~ //int lfil;
//~ double tol = droptol[5];
//~ // BData *w;
//~ //vbiluptr lu = NULL;
//~ vbsptr FF = NULL,FFF = NULL;
//~ vbbsptr Ff,Fff,Bb;
//~ iilutptr luu;
//~
//~ /*---------------------------------------------------------------------
//~ |     Sort the matrix and separate into   |  B  F  |
//~ |                                         |        |
//~ |                                         |  E  C  |
//~ |--------------------------------------------------------------------*/
//~ 
//~ //w = (BData *)Malloc(B->n*sizeof(BData),"main");
//~ //for( i = 0; i < B->n; i++ )
//~ //w[i] = (FLOAT *)Malloc( max_blk_sz, "main" );
//~ //lfil = io.lfil0;
//~ //printf("\n lfil=%d,tol=%20.16e",lfil,tol); getchar();
//~ 
//~ // ierr = vbilutC( B, lu, lfil, tol, w, flog );//int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp )
//~ //ierr = setupVBP4 (vbmat, B->n, C->n, F, E , int *bsz,vbiluptr lu );//int setupVBP4 (vbp4ptr vbmat, int Bn, int Cn,  vbsptr F,  vbsptr E , int *bsz,vbiluptr lu )
//~ FF = (vbsptr) Malloc(sizeof(VBSparMat), "arms2:7" );
//~ Ff = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
//~ Fff = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
//~ Bb = (vbbsptr) Malloc(sizeof(VBBSparMat), "arms2:7" );
//~ luu = (iilutptr) Malloc(sizeof(IIluSpar), "arms2:7" );//luu = (vbbiluptr) Malloc(sizeof(VBBILUSpar), "arms2:7" );
//~ //lu = (vbiluptr)Malloc( sizeof(VBILUSpar), "main" );
//~ /*extern int Bsplit(vbsptr B, vbbsptr Bb, vbsptr F, vbbsptr Ff, int *nBB, int nset);
//~ extern int computLUFinpartition( vbbsptr Bb, vbbsptr Ff, vbbsptr FF, vbbiluptr luu,io_t io, FILE *flog);
//~ extern int Ffassemble( vbsptr FF, vbbsptr Ff);
//~ extern int BLUassemble( vbiluptr lu, vbbiluptr luu);*/
//~ ierr = Bsplit(B, Bb, vbmat->F, Ff, nBB, nset);
//~ ierr = computLUFinpartitionp( Bb, Ff, Fff, luu, droptol, lfil, ipar,flog);
//~ cleanVBBMat(Bb);
//~ cleanVBBMat(Ff);
//~ ierr = Ffassemble( FF, Fff);
//~ cleanVBBMat(Fff);//int cleanVBBMat( vbbsptr vbbmat )
//~ //ierr = BLUassemble( lu, luu);
//~ //cleanVBBILU(luu);//int cleanVBBILU( vbbiluptr luu )
//~ vbmat->plu = luu;
//~
//~ //ierr = computLUF(vbmat->F,vbmat->lu,FF);//int computLUF(vbsptr F, vbiluptr lu, vbsptr FF)
//~ 
//~ //printf("\n     cputimeB^(-1)*F=%f,calendartimeB^(-1)*F=%f\n",tmpp,elapsedtime);
//~ //outputvbmat1(FF,"FFinB.coo",1);
//~ FFF = (vbsptr) Malloc(sizeof(VBSparMat), "arms2:7" );
//~ //outputvbmat1(vbmat->E,"vbmat->E.coo",1);
//~ 
//~ ierr=vbmulvb(1,vbmat->E,FF,FFF);//int vbmulvb(int job, vbsptr A, vbsptr B, vbsptr C)
//~ vbplussvb1(1,C,-1.0,FFF,schur,tol);//int vbplussvb(int job, vbsptr A, FLOAT s, vbsptr B, vbsptr C)??
//~
//~ //outputvbmat(schur,"schur.coo",1);
//~
//~ //for( i = 0; i < B->n; i++ )
//~ //free( w[i] );
//~ //free( w );
//~ cleanVBMat1( FF );
//~ cleanVBMat( FFF );
//~ return 0;
//~ }

int vbscpy(vbsptr amat, vbsptr bmat){
    /*----------------------------------------------------------------------
| Convert CSR matrix to vbSpaFmt struct, copy amat to bmat, vbsptr format, setupbvmat outside??
|----------------------------------------------------------------------
| on entry:
|==========
| ( amat )   = Matrix stored in vbSpaFmt format
|
|
| On return:
|===========
|
| ( bmat )  =  Matrix stored as vbSpaFmt struct containing a copy
|              of amat 
|
|       integer value returned:
|             0   --> successful return.
|             1   --> memory allocation error.
|--------------------------------------------------------------------*/
    int i,j, len, size=amat->n,dimR,dimC,*bsz=amat->bsz,*nB;
    BData *bba;
    int *bja;

    nB = (int *) Malloc(size*sizeof(int), "vbscpy:0" );
    for (j=0; j<size; j++){
        nB[j]=B_DIM(bsz,j);
    }
    if(setupVBMat(bmat,size,nB)) return 1;
    free(nB);
    //printf("\n bmat->bsz=%d,bmat->bszc=%d,bmat->nc=%d,bmat->n=%d",bmat->bsz,bmat->bszc,bmat->nc,bmat->n); getchar();
    /*------------------------------------------------------------*/
    for (i=0; i<size; i++) {
        //printf("\n i=%d",i); getchar();
        dimR = B_DIM(bsz,i);len = bmat->nzcount[i] = amat->nzcount[i];
        if (len > 0) {
            bja = (int *) Malloc(len*sizeof(int), "vbscpy:1" );
            bba = (BData *) Malloc(len*sizeof(BData), "vbscpy:2" );
            for (j=0; j<len; j++){
                //printf("\n j=%d",j); getchar();
                dimC = B_DIM(bsz,amat->ja[i][j]);
                bba[j] = (FLOAT *) Malloc(dimC*dimR*sizeof(FLOAT), "vbscpy:2" );
                //printf("\n bba[j]=%d",bba[j]); getchar();
                copyBData( dimC, dimR, bba[j], amat->ba[i][j], 0 );
                //printf("\n j=%d",j); getchar();
            }
            //printf("\n i=%d",i); getchar();
            memcpy(bja,amat->ja[i],len*sizeof(int));
            //memcpy(bma,amat->ba[i],len*sizeof(FLOAT));
            bmat->ja[i] = bja;
            //printf("\n i=%d",i); getchar();
            bmat->ba[i] = bba;
        }

        else{
            bmat->ja[i] = NULL;
            bmat->ba[i] = NULL;
        }
    }
    return 0;
}

//~ Bug here
//~ int vblusolC( FLOAT *y, FLOAT *x, vbiluptr lu)  
//~ {
//~ /*----------------------------------------------------------------------
//~ *    performs a forward followed by a backward block solve
//~ *    for LU matrix as produced by VBILUT
//~ *    y  = right-hand-side
//~ *    x  = solution on return
//~ *    lu = LU matrix as produced by VBILUT
//~ *
//~ *    note: lu->bf is used to store vector
//~ *--------------------------------------------------------------------*/
//~ int n = lu->n, *bsz = lu->bsz, i, j, icol, dim, sz;
//~ int nzcount, nBs, *ja, inc = 1, OPT;
//~ FLOAT *data, alpha = -1.0, beta = 1.0, alpha2 = 1.0, beta2 = 0.0;
//~ vbsptr L, U;
//~ BData *D, *ba;
//~ 
//~ L = lu->L;
//~ U = lu->U;
//~ D = lu->D;
//~ OPT = lu->DiagOpt;
//~ /* Block L solve */
//~ for( i = 0; i < n; i++ ) {
//~ dim = B_DIM(bsz,i);
//~ nBs = bsz[i];
//~ if (x+nBs != y+nBs)
//~ memcpy(x+nBs, y+nBs, sizeof(FLOAT)*dim);
//~ 
//~ nzcount = L->nzcount[i];
//~ ja = L->ja[i];
//~ ba = L->ba[i];
//~ if (dim == 1)
//~ {
//~ for( j = 0; j < nzcount; j++ ) {
//~ icol = ja[j];
//~ sz = B_DIM(bsz,icol);
//~ if (sz == 1)
//~ x[nBs] -= x[bsz[icol]] * *ba[j];
//~ else
//~ {
//~ data = ba[j];
//~ DGEMV( "n",  dim,  sz,  alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc );
//~ }
//~ }
//~ }
//~ else
//~ {
//~ for( j = 0; j < nzcount; j++ ) {
//~ icol = ja[j];
//~ sz = B_DIM(bsz,icol);
//~ data = ba[j];
//~ DGEMV( "n",  dim,  sz,  alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc );
//~ }
//~ }
//~ }
//~ /* Block -- U solve */
//~ for( i = n-1; i >= 0; i-- ) {
//~ dim = B_DIM(bsz,i);
//~ nzcount = U->nzcount[i];
//~ nBs = bsz[i];
//~ ja = U->ja[i];
//~ ba = U->ba[i];
//~ if (dim == 1)
//~ {
//~ for( j = 0; j < nzcount; j++ ) {
//~ icol = ja[j];
//~ sz = B_DIM(bsz,icol);
//~ if (sz == 1)
//~ x[nBs] -= x[bsz[icol]] * *ba[j];
//~ else
//~ {
//~ data = ba[j];
//~ DGEMV( "n",  dim,  sz,  alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc );
//~ }
//~ }
//~ if (OPT == 1)
//~ x[nBs] /= *D[i];
//~ else
//~ x[nBs] *= *D[i];
//~ }
//~ else
//~ {
//~ for( j = 0; j < nzcount; j++ ) {
//~ icol = ja[j];
//~ sz = B_DIM(bsz,icol);
//~ data = ba[j];
//~ DGEMV( "n", dim, sz, alpha, data, dim, x+bsz[icol], inc, beta, x+nBs, inc );
//~ }
//~ data = D[i];
//~ if (OPT == 1)
//~ luinv( dim, data, x+nBs, lu->bf );
//~ else
//~ DGEMV( "n", dim, dim, alpha2, data, dim, x+nBs, inc, beta2, lu->bf, inc );
//~ memcpy(x+nBs, lu->bf, dim*sizeof(FLOAT));
//~ }
//~ }
//~ 
//~ return 0;
//~ }
