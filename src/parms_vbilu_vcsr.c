#include <math.h>
#include "include/parms_mat_impl.h"
#include "include/parms_opt_impl.h"
#include "DDPQ/globheads.h"
#include "DDPQ/protos.h"

#define DBL_EPSILON 2.2204460492503131e-16 // double epsilon

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

#define qsplit qsplit_ 
#define gauss gauss_
#define zgauss zgauss_
#define bxinv bxinv_
#define zbxinv zbxinv_
//#define dgemm dgemm_

//void dgemm(char*, char*, int*, int*, int*, double*, double*, int*,
//           double*, int*, double*, double*, int*) ;
void gauss (int*, double*, int*);
void zgauss (int*, FLOAT*, int*);
void bxinv (int*, int*, double*, double*, double*);
void zbxinv (int*, int*, FLOAT*, FLOAT*, FLOAT*);
void qsplit(double*, int*, int*, int*) ;

typedef struct parms_vbilu_data {
    vbiluptr lu;
    int schur_start;
    int n;
    int nnz_mat;
    int nnz_prec;
} *parms_vbilu_data;

/*
int qsplitC(FLOAT *a, int *ind, int n, int ncut);
void qqsort(int *ja, FLOAT *ma, int left, int right);
*/

static int parms_vbilu_free(parms_Operator *self)
{
    parms_vbilu_data data;

    data = (parms_vbilu_data)(*self)->data;
    cleanVBILU(data->lu);

    PARMS_FREE(data);
    return 0;
}



static void parms_vbilu_nnz_vcsr(parms_Operator self, int *nnz_mat, int *nnz_pc)
{
    parms_vbilu_data data;

    data = (parms_vbilu_data)self->data;
    *nnz_mat = data->nnz_mat;
    *nnz_pc  = data->nnz_prec;
}

static int parms_vbilu_lsol_vcsr(parms_Operator self, FLOAT *y, FLOAT
                                 *x)
{
    parms_vbilu_data data;
    int start;

    data = (parms_vbilu_data)self->data;
    start = data->schur_start;
    vbLsolp(start, data->lu, y, x);

    return 0;
}

static int parms_vbilu_invs_vcsr(parms_Operator self,  FLOAT *y, FLOAT
                                 *x)
{
    parms_vbilu_data data;
    int start;

    data = (parms_vbilu_data)self->data;
    start = data->schur_start;
    vbinvsp(start, data->lu, y, x);
    return 0;
}

static int parms_vbilu_ascend_vcsr(parms_Operator self, FLOAT *y, FLOAT
                                   *x)
{
    parms_vbilu_data data;
    int start;

    data = (parms_vbilu_data)self->data;
    start = data->schur_start;
    vbUsolp(start, data->lu, y, x);
    return 0;
}

static int parms_vbilu_sol_vcsr(parms_Operator self, FLOAT *y,
                                FLOAT *x)
{
    parms_vbilu_data data;

    data = (parms_vbilu_data)self->data;
    vblusolC(y, x, data->lu);
    return 0;
}

static int parms_vbilu_getssize_vcsr(parms_Operator self)
{
    parms_vbilu_data data;

    data = (parms_vbilu_data)self->data;
    return data->schur_start;
}

static struct parms_Operator_ops parms_vbilu_sol_vptr = {
    parms_vbilu_sol_vcsr,
    parms_vbilu_lsol_vcsr,
    parms_vbilu_invs_vcsr,
    parms_vbilu_ascend_vcsr,
    parms_vbilu_getssize_vcsr,
    parms_vbilu_nnz_vcsr,
    parms_vbilu_free,
    NULL//parms_vbilu_view
};

int parms_vbilut_vcsr(parms_Mat self, parms_FactParam param, void *mat,
                      parms_Operator *op)
{
    parms_Operator newOpt;
    parms_vbilu_data data;
    parms_bvcsr     mat_bvcsr;
    parms_Map      is;
    int i, j, j1, j2, n, m, start, schur_start, max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE;
    int lenl, lenu, len, k, jpos, jrow, ierr;
    int *bsz;
    FLOAT one = 1.0, zero = 0.0;
    int dim, szjrow, sz, col, upos, para, kk, lfil;
    int nzcount, *ja, *jbuf, *iw;
    double t, tnorm, tolnorm, *xnrm, *wn, vbnorm2( int, FLOAT * ), tol;
    vbsptr L, U;
    BData *ba, *D, buf_fact, buf_ns;

//#if defined(DBL_CMPLX)
//    double shf = 0.0, ti, sgny;
//    int nnzA = 0;
//    /*----------get nnz of mat ---------*/
//    for(j = 0; j<((parms_bvcsr)mat)->n; j++)
//        nnzA += ((parms_bvcsr)mat)->nzcount[j];
//#endif

    is = self->is;
    n = param->n;
    start = param->start;
    schur_start = param->schur_start;
    printf("%d %d\n", start, schur_start);

    if (schur_start == -1) {
        schur_start = n;
    }

    mat_bvcsr = (parms_bvcsr)mat;
    printf("mat_bvcsr->bsz[is->schur_start]=%d\n", mat_bvcsr->bsz[is->schur_start]);

    //exit(1);

    if (!param->isalloc) {
        parms_OperatorCreate(&newOpt);
        PARMS_MEMCPY(newOpt->ops, &parms_vbilu_sol_vptr, 1);
        PARMS_NEW(data);
        PARMS_NEW(data->lu);
        setupVBILU(data->lu, n, mat_bvcsr->bsz);
        data->nnz_mat  = 0;
        data->nnz_prec = 0;
        param->isalloc = true;

        newOpt->data = data;
        *op = newOpt;
    }
    else {
        data = (*op)->data;
    }
    BData *w;
    PARMS_NEWARRAY(w, n);
    for( i = 0; i < n; i++ )
        PARMS_NEWARRAY(w[i], max_blk_sz);

    L = data->lu->L;
    U = data->lu->U;
    D = data->lu->D;
    bsz = mat_bvcsr->bsz;

    //~ vbilutC(mat_bvcsr, data->lu, param->lfil[0], param->droptol[0], w, stdout);

    /* malloc memory for working arrays */
    iw = data->lu->work;
    PARMS_NEWARRAY(jbuf, n);
    PARMS_NEWARRAY(buf_fact, max_blk_sz);
    PARMS_NEWARRAY(buf_ns, max_blk_sz);
    PARMS_NEWARRAY(xnrm, n);
    PARMS_NEWARRAY(wn, n);

    //~ goto done;

    /* set indicator array jw to -1 */
    for( i = 0; i < n; i++ ) iw[i] = -1;

    tol = param->droptol[0];
    lfil = param->lfil[0];

    /* beginning of main loop */
    ierr = 0;
    m = mat_bvcsr->n;
    for( i = 0; i < m; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        nzcount = mat_bvcsr->nzcount[i];
        ja = mat_bvcsr->ja[i];
        ba = mat_bvcsr->ba[i];
        j1 = 0;
        j2 = nzcount;
        tnorm = 0.0;
        for( j = 0; j < nzcount; j++ ) {
            sz = B_DIM(bsz,ja[j]);
            t = vbnorm2( dim * sz, ba[j] );
            tnorm = max(t,tnorm);
        }
        if( tnorm == 0.0 ) {
            printf("vbilut:  zero row encountered.\n" );
            // go to done
            return -2;
        }
        tolnorm = tol * tnorm;
        /* unpack L-part and U-part of row of A in array w */
        lenu = 0;
        lenl = 0;
        if (i+start < schur_start) {
            jbuf[i+start] = i+start;
            zrmC( dim, dim, w[i+start] );
            iw[i+start] = i+start;
            lenu = 1;
            for (j = j1; j < j2; j++) {
                col = ja[j];
                sz = B_DIM(bsz,col);
                t = vbnorm2( dim * sz, ba[j] );
                if( t < tolnorm && col != i + start ) continue;
                if( col < i + start) {
                    jbuf[lenl] = col;
                    iw[col] = lenl;
                    copyBData( dim, sz, w[lenl], ba[j], 0 );
                    lenl++;
                }
                else if(col == i+start){
                    copyBData( dim, dim, w[i+start], ba[j], 0 );
                }
                else {
                    jpos = i + lenu + start;
                    jbuf[jpos] = col;
                    iw[col] = jpos;
                    copyBData( dim, sz, w[jpos], ba[j], 0 );
                    lenu++;
                }
            }
        }
        else {
            for (j = j1; j < j2; j++) {
                col = ja[j];
                sz = B_DIM(bsz,col);
                t = vbnorm2( dim * sz, ba[j] );
                if( t < tolnorm && col != schur_start ) continue;
                if( col < schur_start) {
                    jbuf[lenl] = col;
                    iw[col] = lenl;
                    copyBData( dim, sz, w[lenl], ba[j], 0 );
                    lenl++;
                }
                else {
                    jpos = lenu + schur_start;
                    jbuf[jpos] = col;
                    iw[col] = jpos;
                    copyBData( dim, sz, w[jpos], ba[j], 0 );
                    lenu++;
                }
            }
        }

        j = -1;
        len = 0;
        /* eliminate previous rows */
        while( ++j < lenl ) {
            /*----------------------------------------------------------------------
 *   in order to do the elimination in the correct order we must select
 *   the smallest column index among jw(k), k=j+1, ..., lenl.
 *--------------------------------------------------------------------*/
            jrow = jbuf[j];
            jpos = j;
            /* determine smallest column index */
            for( k = j + 1; k < lenl; k++ ) {
                if( jbuf[k] < jrow ) {
                    jrow = jbuf[k];
                    jpos = k;
                }
            }
            szjrow = B_DIM(bsz,jrow);
            if( jpos != j ) {
                col = jbuf[j];
                jbuf[j] = jbuf[jpos];
                jbuf[jpos] = col;
                iw[jrow] = j;
                iw[col] = jpos;
                sz = B_DIM(bsz,col);
                copyBData( dim, sz, buf_ns, w[j], 0 );
                copyBData( dim, szjrow, w[j], w[jpos], 0 );
                copyBData( dim, sz, w[jpos], buf_ns, 0 );
            }
            /* get the multiplier for row to be eliminated (jrow). */
            /* fact = w(jj)*alu(jrow)                              */
#if defined(DBL_CMPLX)
            zbxinv( &dim, &szjrow, D[jrow], w[j], buf_fact );
#else
            bxinv( &dim, &szjrow, D[jrow], w[j], buf_fact );
#endif
            /* zero out element in row by resetting iw(jrow) to -1 */
            iw[jrow] = -1;

            if( vbnorm2( dim * szjrow, buf_fact ) * xnrm[jrow] <= tolnorm )
                continue;

            /* combine current row and row jrow */
            nzcount = U->nzcount[jrow];
            ja = U->ja[jrow];
            ba = U->ba[jrow];
            if (i+start < schur_start) { /* ilut */
                for( k = 0; k < nzcount; k++ ) {
                    col = ja[k];
                    sz = B_DIM(bsz,col);
                    //      dgemm ("n", "n", &dim, &sz, &szjrow, &one, buf_fact, &dim,
                    //		 ba[k], &szjrow, &zero, buf_ns, &dim );
                    GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                           ba[k], szjrow, zero, buf_ns, dim );

                    jpos = iw[col];

                    /* if fill-in element is small then disregard: */
                    if( vbnorm2( dim * sz, buf_ns ) < tolnorm && jpos == -1 ) continue;

                    if( col >= i + start) {
                        /* dealing with upper part */
                        if( jpos == -1 && vbnorm2( dim * sz, buf_ns ) > tolnorm) {
                            /* this is a fill-in element */
                            upos = i + lenu + start;
                            jbuf[upos] = col;
                            iw[col] = upos;
                            copyBData( dim, sz, w[upos], buf_ns, 0 );
                            ++lenu;
                        }
                        else {
                            /* this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                        }
                    }
                    else {
                        /* dealing with lower part */
                        if(jpos == -1) {
                            /* this is a fill-in element */
                            jbuf[lenl] = col;
                            iw[col] = lenl;
                            copyBData( dim, sz, w[lenl], buf_ns, 0 );
                            lenl++;
                        } else {
                            /*-------------------- this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                        }
                    }
                }
            }
            else {		/* partial ilut */
                for (k = 0; k < nzcount; k++) {
                    col = ja[k];
                    sz = B_DIM(bsz,col);
                    //	  dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                    //		 ba[k], &szjrow, &zero, buf_ns, &dim );
                    GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                           ba[k], szjrow, zero, buf_ns, dim );
                    jpos = iw[col];

                    /* if fill-in element is small then disregard: */
                    if( vbnorm2( dim * sz, buf_ns ) < tolnorm && jpos == -1 ) continue;

                    if( col >= schur_start) {
                        /* dealing with upper part */
                        if( jpos == -1 && vbnorm2( dim * sz, buf_ns ) > tolnorm) {
                            /* this is a fill-in element */
                            upos = schur_start + lenu;
                            jbuf[upos] = col;
                            iw[col] = upos;
                            copyBData( dim, sz, w[upos], buf_ns, 0 );
                            ++lenu;
                        }
                        else {
                            /* this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                        }
                    }
                    else {
                        /* dealing with lower part */
                        if(jpos == -1) {
                            /* this is a fill-in element */
                            jbuf[lenl] = col;
                            iw[col] = lenl;
                            copyBData( dim, sz, w[lenl], buf_ns, 0 );
                            lenl++;
                        } else {
                            /*-------------------- this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                        }
                    }
                }
            }

            /* store this pivot element -- (from left to right -- no danger
      of overlap with the working elements in L (pivots). */
            copyBData( dim, szjrow, w[len], buf_fact, 1 );
            jbuf[len] = jrow;
            len++;
        }

        /* update l-matrix */
        lenl = len;
        len = min( lenl, lfil );
        for( j = 0; j < lenl; j++ ) {
            sz = B_DIM(bsz,jbuf[j]);
            wn[j] = vbnorm2( dim*sz, w[j] );
            iw[j] = j;
        }
        qsplit( wn, iw, &lenl, &len );
        L->nzcount[i+start] = len;
        if( len > 0 ) {
            PARMS_NEWARRAY(L->ja[i+start], len);
            PARMS_NEWARRAY(L->ba[i+start], len);
        }
        ja = L->ja[i+start];
        ba = L->ba[i+start];
        for( j = 0; j < len; j++ ) {
            jpos = iw[j];
            ja[j] = jbuf[jpos];
            sz = B_DIM(bsz,ja[j]);
            PARMS_NEWARRAY(ba[j], dim*sz);
            copyBData( dim, sz, ba[j], w[jpos], 0 );
        }
        for( j = 0; j < lenl; j++ ) iw[j] = -1;

        /* update u-matrix */
        if (i+start < schur_start) {
            len = min( lenu, lfil );
            for( j = 1; j < lenu; j++ ) {
                jpos = i+start+j;
                sz = B_DIM(bsz,jbuf[jpos]);
                wn[j-1] = vbnorm2( dim*sz, w[jpos] );
                iw[j-1] = jpos;
            }
            para = lenu - 1;
            qsplit( wn, iw, &para, &len );
            nzcount = U->nzcount[i+start] = len-1;
            if( nzcount > 0 ) {
                PARMS_NEWARRAY(U->ja[i+start], nzcount);
                PARMS_NEWARRAY(U->ba[i+start], nzcount);
            }
            ja = U->ja[i+start];
            ba = U->ba[i+start];
            t = vbnorm2( dim*dim, w[i+start] );
            for( j = 0; j < nzcount; j++ ) {
                jpos = iw[j];
                ja[j] = jbuf[jpos];
                sz = B_DIM(bsz,ja[j]);
                PARMS_NEWARRAY(ba[j], dim*sz);
                copyBData( dim, sz, ba[j], w[jpos], 0 );
                t = max(t,wn[j]);
            }
            xnrm[i+start] = t;

            /* store inverse of diagonal element of u */
            PARMS_NEWARRAY(D[i+start], dim*dim);
            copyBData( dim, dim, D[i+start], w[i+start], 0 );

#if defined(DBL_CMPLX)
            zgauss( &dim, D[i+start], &ierr );
#else
            gauss( &dim, D[i+start], &ierr );
#endif
        }
        else if (i >= schur_start) {
            U->nzcount[start+i] = lenu-1;
            PARMS_NEWARRAY(U->ja[start+i], lenu);
            PARMS_NEWARRAY(U->ba[start+i], lenu);
            ja = U->ja[i+start];
            ba = U->ba[i+start];
            for( j = 0; j < lenu; j++ ) {
                jpos = iw[j];
                ja[j] = jbuf[jpos];
                sz = B_DIM(bsz,ja[j]);
                PARMS_NEWARRAY(ba[j], dim*sz);
                copyBData( dim, sz, ba[j], w[jpos], 0 );
            }
        }
        for( j = 0; j < lenu; j++ ) iw[j] = -1;
        if (i+start < schur_start)
            for( j = 0; j < lenu; j++ ) iw[jbuf[i+start+j]] = -1;
        else
            for( j = 0; j < lenu; j++ ) iw[jbuf[schur_start+j]] = -1;
    }

    //done:

    data->lu->DiagOpt = 1;
    printf("is->schur_start=%d\n", is->schur_start);
    //printf("bsz[is->schur_start]=%d\n", bsz[is->schur_start]);
    //output_intvectorpa("bszin_vbilut",bsz,0, n+1);

    //exit(1);

    data->schur_start = bsz[is->schur_start];
    data->n = bsz[n];

    for( i = 0; i < n; i++ )
        PARMS_FREE(w[i]);
    PARMS_FREE(w);
    PARMS_FREE(jbuf);
    PARMS_FREE(buf_fact);
    PARMS_FREE(buf_ns);
    PARMS_FREE(xnrm);
    PARMS_FREE(wn);

    /* compute the number of nonzeros in matrix */
    data->nnz_mat = nnzVBMat1(mat_bvcsr);

    /* computer the number of nonzeros in pc */
    data->nnz_prec = nnz_vbilu(data->lu);

    return ierr;
}

int parms_vbiluk_vcsr(parms_Mat self, parms_FactParam param, void *mat,
                      parms_Operator *op)
{
    parms_Operator newOpt;
    parms_vbilu_data data;
    parms_bvcsr     mat_bvcsr;
    parms_Map      is;
    int i, j, n, m, start, schur_start, max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE;
    int lenl, lenu, len, k, jpos, jrow, ierr;
    int *bsz;
    FLOAT one = 1.0, zero = 0.0;
    int dim, szjrow, sz, col, upos, kk, lfil;
    int nzcount, *ja, *jbuf, *iw, *levls, **levs, *lev, jlev;
    vbsptr L, U;
    BData *ba, *D, buf_fact, buf_ns;

    //~ #if defined(DBL_CMPLX)
    //~ double shf = 0.0, ti, sgny;
    //~ int nnzA = 0;
    //~ /*----------get nnz of mat ---------*/
    //~ for(j = 0; j<((parms_bvcsr)mat)->n; j++)
    //~ nnzA += ((parms_bvcsr)mat)->nnzrow[j];
    //~ #endif

    is = self->is;
    n = param->n;
    start = param->start;
    schur_start = param->schur_start;
    printf("%d %d\n", start, schur_start);

    if (schur_start == -1) {
        schur_start = n;
    }

    mat_bvcsr = (parms_bvcsr)mat;

    if (!param->isalloc) {
        parms_OperatorCreate(&newOpt);
        PARMS_MEMCPY(newOpt->ops, &parms_vbilu_sol_vptr, 1);
        PARMS_NEW(data);
        PARMS_NEW(data->lu);
        setupVBILU(data->lu, n, mat_bvcsr->bsz);
        data->nnz_mat  = 0;
        data->nnz_prec = 0;
        param->isalloc = true;

        newOpt->data = data;
        *op = newOpt;
    }
    else {
        data = (*op)->data;
    }
    BData *w;
    PARMS_NEWARRAY(w, n);
    for( i = 0; i < n; i++ )
        PARMS_NEWARRAY(w[i], max_blk_sz);

    L = data->lu->L;
    U = data->lu->U;
    D = data->lu->D;
    bsz = mat_bvcsr->bsz;

    //~ vbilukC(param->lfil[0], mat_bvcsr, data->lu, stdout);
    //~ goto done;

    /* malloc memory for working arrays */
    iw = data->lu->work;
    PARMS_NEWARRAY(jbuf, n);
    PARMS_NEWARRAY(buf_fact, max_blk_sz);
    PARMS_NEWARRAY(buf_ns, max_blk_sz);
    PARMS_NEWARRAY(levls, n);
    PARMS_NEWARRAY(levs, n);

    /* set indicator array jw to -1 */
    for( i = 0; i < n; i++ ) {
        iw[i] = -1;
        levls[i] = 0;
    }

    //lfil = param->lfil[0];
    lfil = param->ipar[0];

    /* beginning of main loop */
    ierr = 0;
    m = mat_bvcsr->n;
    for( i = 0; i < m; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        nzcount = mat_bvcsr->nzcount[i];
        ja = mat_bvcsr->ja[i];
        ba = mat_bvcsr->ba[i];
        //~ tnorm = 0.0;
        //~ for( j = 0; j < nzcount; j++ ) {
        //~ sz = B_DIM(bsz,ja[j]);
        //~ t = vbnorm2( dim * sz, ba[j] );
        //~ tnorm = max(t,tnorm);
        //~ }
        //~ if( tnorm == 0.0 ) {
        //~ printf("vbilut:  zero row encountered.\n" );
        //~ // go to done
        //~ return -2;
        //~ }
        /* unpack L-part and U-part of row of A in array w */
        lenu = 0;
        lenl = 0;
        if (i+start < schur_start) {
            //TODO: Add a shift thing (see pARMS iluk)
            jbuf[i+start] = i+start;
            zrmC( dim, dim, w[i+start] );
            iw[i+start] = i+start;
            lenu = 1;
            for (j = 0; j < nzcount; j++) {
                col = ja[j];
                sz = B_DIM(bsz,col);
                if( col < i + start) {
                    jbuf[lenl] = col;
                    levls[lenl] = 0;
                    iw[col] = lenl;
                    copyBData( dim, sz, w[lenl], ba[j], 0 );
                    lenl++;
                }
                else if(col == i+start){
                    levls[i+start] = 0;
                    copyBData( dim, dim, w[i+start], ba[j], 0 );
                }
                else {
                    jpos = i + lenu + start;
                    jbuf[jpos] = col;
                    levls[jpos] = 0;
                    iw[col] = jpos;
                    copyBData( dim, sz, w[jpos], ba[j], 0 );
                    lenu++;
                }
            }
        }
        else {
            for (j = 0; j < nzcount; j++) {
                col = ja[j];
                sz = B_DIM(bsz,col);
                if( col < schur_start) {
                    jbuf[lenl] = col;
                    levls[lenl] = 0;
                    iw[col] = lenl;
                    copyBData( dim, sz, w[lenl], ba[j], 0 );
                    lenl++;
                }
                else {
                    jpos = lenu + schur_start;
                    jbuf[jpos] = col;
                    levls[jpos] = 0;
                    iw[col] = jpos;
                    copyBData( dim, sz, w[jpos], ba[j], 0 );
                    lenu++;
                }
            }
        }

        j = -1;
        /* eliminate previous rows */
        while( ++j < lenl ) {
            /*----------------------------------------------------------------------
 *   in order to do the elimination in the correct order we must select
 *   the smallest column index among jw(k), k=j+1, ..., lenl.
 *--------------------------------------------------------------------*/
            jrow = jbuf[j];
            jpos = j;
            /* determine smallest column index */
            for( k = j + 1; k < lenl; k++ ) {
                if( jbuf[k] < jrow ) {
                    jrow = jbuf[k];
                    jpos = k;
                }
            }
            szjrow = B_DIM(bsz,jrow);
            if( jpos != j ) {
                col = jbuf[j];
                jbuf[j] = jbuf[jpos];
                jbuf[jpos] = col;
                iw[jrow] = j;
                iw[col] = jpos;
                sz = B_DIM(bsz,col);
                col = levls[j];
                levls[j] = levls[jpos];
                levls[jpos] = col;
                copyBData( dim, sz, buf_ns, w[j], 0 );
                copyBData( dim, szjrow, w[j], w[jpos], 0 );
                copyBData( dim, sz, w[jpos], buf_ns, 0 );
            }
            /* get the multiplier for row to be eliminated (jrow). */
            /* fact = w(jj)*alu(jrow)                              */
#if defined(DBL_CMPLX)
            zbxinv( &dim, &szjrow, D[jrow], w[j], buf_fact );
#else
            bxinv( &dim, &szjrow, D[jrow], w[j], buf_fact );
#endif
            /* zero out element in row by resetting iw(jrow) to -1 */
            iw[jrow] = -1;

            /* combine current row and row jrow */
            nzcount = U->nzcount[jrow];
            ja = U->ja[jrow];
            ba = U->ba[jrow];
            lev = levs[jrow];
            jlev = levls[j];

            if (jlev > lfil) continue;

            if (i+start < schur_start) {
                for( k = 0; k < nzcount; k++ ) {
                    col = ja[k];
                    sz = B_DIM(bsz,col);
                    //	  dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                    //		 ba[k], &szjrow, &zero, buf_ns, &dim );
                    GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                           ba[k], szjrow, zero, buf_ns, dim );
                    jpos = iw[col];

                    if( col >= i + start) {
                        /* dealing with upper part */
                        if( jpos == -1) {
                            /* this is a fill-in element */
                            upos = i + lenu + start;
                            jbuf[upos] = col;
                            levls[upos] = jlev+lev[k]+1;
                            iw[col] = upos;
                            copyBData( dim, sz, w[upos], buf_ns, 0 );
                            ++lenu;
                        }
                        else {
                            /* this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                            levls[jpos] = min(levls[jpos], jlev+lev[k]+1);
                        }
                    }
                    else {
                        /* dealing with lower part */
                        if(jpos == -1) {
                            /* this is a fill-in element */
                            jbuf[lenl] = col;
                            levls[lenl] = jlev+lev[k]+1;
                            iw[col] = lenl;
                            copyBData( dim, sz, w[lenl], buf_ns, 0 );
                            lenl++;
                        } else {
                            /*-------------------- this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                            levls[jpos] = min(levls[jpos], jlev+lev[k]+1);
                        }
                    }
                }
            }
            else {
                for (k = 0; k < nzcount; k++) {
                    col = ja[k];
                    sz = B_DIM(bsz,col);
                    GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                           ba[k], szjrow, zero, buf_ns, dim );
                    //      dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                    //             ba[k], &szjrow, &zero, buf_ns, &dim );
                    jpos = iw[col];

                    if( col >= schur_start) {
                        /* dealing with upper part */
                        if( jpos == -1) {
                            /* this is a fill-in element */
                            upos = schur_start + lenu;
                            jbuf[upos] = col;
                            levls[upos] = jlev+lev[k]+1;
                            iw[col] = upos;
                            copyBData( dim, sz, w[upos], buf_ns, 0 );
                            ++lenu;
                        }
                        else {
                            /* this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                            levls[jpos] = min(levls[jpos], jlev+lev[k]+1);
                        }
                    }
                    else {
                        /* dealing with lower part */
                        if(jpos == -1) {
                            /* this is a fill-in element */
                            jbuf[lenl] = col;
                            levls[lenl] = jlev+lev[k]+1;
                            iw[col] = lenl;
                            copyBData( dim, sz, w[lenl], buf_ns, 0 );
                            lenl++;
                        } else {
                            /*-------------------- this is not a fill-in element */
                            for (kk=0;kk<dim*sz; kk++)
                                w[jpos][kk] +=  buf_ns[kk];
                            levls[jpos] = min(levls[jpos], jlev+lev[k]+1);
                        }
                    }
                }
            }

            /* store this pivot element -- (from left to right -- no danger
      of overlap with the working elements in L (pivots). */
            copyBData( dim, szjrow, w[j], buf_fact, 1 );
            jbuf[j] = jrow;
        }

        /* update l-matrix */
        len = 0;
        for( j = 0; j < lenl; j++ ) {
            if (levls[j] <= lfil) {
                iw[len] = j;
                len++;
            }
        }
        L->nzcount[i+start] = len;
        if( len > 0 ) {
            PARMS_NEWARRAY(L->ja[i+start], len);
            PARMS_NEWARRAY(L->ba[i+start], len);
        }
        ja = L->ja[i+start];
        ba = L->ba[i+start];
        for( j = 0; j < len; j++ ) {
            jpos = iw[j];
            ja[j] = jbuf[jpos];
            sz = B_DIM(bsz,ja[j]);
            PARMS_NEWARRAY(ba[j], dim*sz);
            copyBData( dim, sz, ba[j], w[jpos], 0 );
        }
        for( j = 0; j < lenl; j++ ) iw[j] = -1;

        /* update u-matrix */
        if (i+start < schur_start) {
            len = 0;
            for( j = 1; j < lenu; j++ ) {
                jpos = i+start+j;
                if (levls[jpos] <= lfil)
                {
                    iw[len] = jpos;
                    len++;
                }
            }
        }
        else
        {
            len = 0;
            for( j = 0; j < lenu; j++ ) {
                jpos = schur_start+j;
                if (levls[jpos] <= lfil)
                {
                    iw[len] = jpos;
                    len++;
                }
            }
        }

        nzcount = U->nzcount[i+start] = len;
        if( nzcount > 0 ) {
            PARMS_NEWARRAY(U->ja[i+start], nzcount);
            PARMS_NEWARRAY(U->ba[i+start], nzcount);
            PARMS_NEWARRAY(levs[i+start], nzcount);
        }
        ja = U->ja[i+start];
        ba = U->ba[i+start];
        lev = levs[i+start];

        for( j = 0; j < nzcount; j++ ) {
            jpos = iw[j];
            ja[j] = jbuf[jpos];
            lev[j] = levls[jpos];
            sz = B_DIM(bsz,ja[j]);
            PARMS_NEWARRAY(ba[j], dim*sz);
            copyBData( dim, sz, ba[j], w[jpos], 0 );
        }

        if (i+start < schur_start)
        {
            /* store inverse of diagonal element of u */
            PARMS_NEWARRAY(D[i+start], dim*dim);
            copyBData( dim, dim, D[i+start], w[i+start], 0 );

#if defined(DBL_CMPLX)
            zgauss( &dim, D[i+start], &ierr );
#else
            gauss( &dim, D[i+start], &ierr );
#endif
        }
        for( j = 0; j < lenu; j++ ) iw[j] = -1;
        if (i+start < schur_start)
            for( j = 0; j < lenu; j++ ) iw[jbuf[i+start+j]] = -1;
        else
            for( j = 0; j < lenu; j++ ) iw[jbuf[schur_start+j]] = -1;
    }

    data->lu->DiagOpt = 1;
    data->schur_start = bsz[is->schur_start];
    data->n = bsz[n];

    for (i=0; i < m-1; i++)
        if (U->nzcount[i+start])
            PARMS_FREE(levs[i+start]);
    PARMS_FREE(levs);

    for( i = 0; i < n; i++ )
        PARMS_FREE(w[i]);
    PARMS_FREE(w);
    PARMS_FREE(jbuf);
    PARMS_FREE(buf_fact);
    PARMS_FREE(buf_ns);
    PARMS_FREE(levls);

    //~ done:

    /* compute the number of nonzeros in matrix */
    data->nnz_mat = nnzVBMat1(mat_bvcsr);

    /* computer the number of nonzeros in pc */
    data->nnz_prec = nnz_vbilu(data->lu);

    return ierr;
}
