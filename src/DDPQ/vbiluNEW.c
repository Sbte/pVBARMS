#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "globheads.h"
#include "protos.h"

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

void bxinv (int*, int*, double*,double*,double*);
void zbxinv (int*, int*, FLOAT*,FLOAT*,FLOAT*);
void qsplit(double*, int*, int*, int*) ;
void gauss (int *, double*, int*);
void zgauss (int *, FLOAT*, int*);

int vbiluNEW(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol, 
             int *lfil, vbsptr schur)
{
    /*----------------------------------------------------------------------------
 * Block ILUT (BILUT) preconditioner
 * Block incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal blocks
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * vbmat    = block matrix stored in VBSpaFmt format -- see globheads.h for
 *            details on format, the block sizes might be different
 * lu       = pointer to a VBILUSpar struct -- see globheads.h for details
 *            on format
 * lfil     = integer. The fill-in parameter. Each row of L and
 *            each row of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS.
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * w        = working array
 * fp       = file pointer for error log ( might be stdout )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> singular block or zero row encountered
 * lu->n    = dimension of the block matrix
 *   ->bsz  = the row/col of the first element of each diagonal block
 *            the size of the i-th row block should be bsz[i+1] - bsz[i]
 *   ->L    = L part -- stored in VBSpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in VBSpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonal blocks of the input block matrix must not be singular
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows.
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th row of L
 *    and the largest lfil elements in the i-th row of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/
    int n = schur->n, *bsz = vbmat->bsz, *bszc = vbmat->F->bszc, ierr;
    FLOAT one = 1.0, zero = 0.0;
    int dim, szjrow, sz, len, lenu, lenl, len2, col, jpos, jrow, upos, para;
    int max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(double);
    int lnzcount, rnzcount, *jbuf, *iw, i, j, k, kk, *lja, *rja;
    int *jbuf2, *iw2, lsize, rsize, **lfja, *lflen;
    int fil0=lfil[0],fil1=lfil[1],fil2=lfil[2],fil4=lfil[4];
    double tnorm;
    double drop0=droptol[0], drop1=droptol[1], drop2=droptol[2];
    double drop3=droptol[3], drop4=droptol[4];
    double t, *xnrm, *wn, vbnorm2( int, FLOAT * );
    vbsptr L, U;
    BData *lba, *rba, *D, buf_fact, buf_ns;
    BData *w, *w2, **lfma;

    lsize = vbmat->nB;
    rsize = C->n;
    printf("lsize %d, rsize %d, n %d\n", lsize, rsize, n);
    n = lsize > rsize ? lsize : rsize;

    if(fil0 < 0 || fil1 < 0 || n <= 0) {
        printf( "vbiluNEW: Illegal value for lfil.\n" );
        return -1;
    }
    buf_fact = (BData)Malloc( max_blk_sz, "vbiluNEW:1" );
    buf_ns = (BData)Malloc( max_blk_sz, "vbiluNEW:2" );
    xnrm = (double *)Malloc( n * sizeof(double), "vbiluNEW:3" );
    wn = (double *)Malloc( n * sizeof(double), "vbiluNEW:4" );

    jbuf = (int *)Malloc( n * sizeof(int), "vbiluNEW:5" );
    w = (BData *)Malloc( n * sizeof(BData),"vbiluNEW:6");
    for( i = 0; i < n; i++ )
        w[i] = (FLOAT *)Malloc( max_blk_sz, "vbiluNEW:7" );
    iw = (int *)Malloc( n * sizeof(int), "vbiluNEW:8" );
    //~ lu = (vbiluptr)Malloc( sizeof(VBILUSpar), "vbiluNEW" );
    //~ lu->DiagOpt = 0;
    //~
    //~ setupVBILU( lu, n, bsz );
    L = vbmat->lu->L;
    U = vbmat->lu->U;
    D = vbmat->lu->D;
    //~
    //~ iw = lu->work;
    //~
    jbuf2 = (int *)Malloc( n * sizeof(int), "vbiluNEW:9" );
    w2 = (BData *)Malloc( n * sizeof(BData),"vbiluNEW:10");
    for( i = 0; i < n; i++ )
        w2[i] = (FLOAT *)Malloc( max_blk_sz, "vbiluNEW:11" );
    iw2 = (int *)Malloc( n * sizeof(int), "vbiluNEW:12" );

    lfma = (BData **) Malloc(lsize*sizeof(BData *), "vbiluNEW:13" );
    lfja = (int **) Malloc(lsize*sizeof(int *), "vbiluNEW:14" );
    lflen = (int *) Malloc(lsize*sizeof(int), "vbiluNEW:15" );
    //~ lu2 = (vbiluptr)Malloc( sizeof(VBILUSpar), "vbiluNEW" );
    //~ lu2->DiagOpt = 0;
    //~
    //~ setupVBILU( lu2, n, bsz );
    //~ L2 = lu->L;
    //~ U2 = lu->U;
    //~ D2 = lu->D;
    //~
    //~ iw2 = lu->work;

    /* set indicator array jw to -1 */
    for( i = 0; i < n; i++ )
    {
        iw[i] = -1;
        iw2[i] = -1;
    }

    /* beginning of main loop */
    for( i = 0; i < lsize; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        lnzcount = B->nzcount[i];
        lja = B->ja[i];
        lba = B->ba[i];
        rnzcount = vbmat->F->nzcount[i];
        rja = vbmat->F->ja[i];
        rba = vbmat->F->ba[i];
        tnorm = 0.0;
        for( j = 0; j < lnzcount; j++ ) {
            sz = B_DIM(bsz,lja[j]);
            t = vbnorm2( dim * sz, lba[j] );
            tnorm = max(t,tnorm);
        }
        if( tnorm == 0.0 ) {
            printf( "vbiluNEW:  zero row encountered.\n" );
            return -2;
        }
        //~ tnorm = 0.0;
        // Something complex ???

        /* unpack L-part and U-part of row of A in arrays w */
        lenu = 1;
        lenl = 0;
        jbuf[i] = i;
        zrmC( dim, dim, w[i] );
        iw[i] = i;

        for( j = 0; j < lnzcount; j++ ) {
            col = lja[j];
            sz = B_DIM(bsz, col);
            t = vbnorm2( dim * sz, lba[j] );
            if( t < tnorm * drop0 && col != i ) continue;
            if( col < i ) {
                jbuf[lenl] = col;
                iw[col] = lenl;
                copyBData( dim, sz, w[lenl], lba[j], 0 );
                lenl++;
            } else if( col == i ) {
                copyBData( dim, dim, w[i], lba[j], 0 );
            } else {
                jpos = i + lenu;
                jbuf[jpos] = col;
                iw[col] = jpos;
                copyBData( dim, sz, w[jpos], lba[j], 0 );
                lenu++;
            }
        }

        len2 = 0;
        for ( j = 0; j < rnzcount; j++ ) {
            col = rja[j];
            sz = B_DIM(bszc, col);
            jbuf2[len2] = col;
            //~ printf("%d %d %d %d\n", i, j, n, col);
            copyBData( dim, sz, w2[len2], rba[j], 0 );
            iw2[col] = len2;
            len2++;
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
#endif                /* zero out element in row by resetting jw(n+jrow) to -1 */
            iw[jrow] = -1;

            if( vbnorm2( dim * szjrow, buf_fact ) * xnrm[jrow] <= tnorm * drop0 )
                continue;

            /* combine current row and row jrow */
            lnzcount = U->nzcount[jrow];
            lja = U->ja[jrow];
            lba = U->ba[jrow];
            rnzcount = lflen[jrow];
            rja = lfja[jrow];
            rba = lfma[jrow];
            for( k = 0; k < lnzcount; k++ ) {
                col = lja[k];
                sz = B_DIM(bsz,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       lba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       lba[k], szjrow, zero, buf_ns, dim );
                jpos = iw[col];

                /* if fill-in element is small then disregard: */
                if( vbnorm2( dim * sz, buf_ns ) < tnorm * drop1 && jpos == -1 ) continue;

                if( col >= i ) {
                    /* dealing with upper part */
                    //          if( jpos == -1 ) {
                    if( jpos == -1 && vbnorm2( dim * sz, buf_ns ) > tnorm * drop1) {
                        /* this is a fill-in element */
                        upos = i + lenu;
                        jbuf[upos] = col;
                        iw[col] = upos;
                        copyBData( dim, sz, w[upos], buf_ns, 0 );
                        lenu++;
                    } else {
                        /*-------------------- this is not a fill-in element */
                        for (kk=0;kk<dim*sz; kk++)
                            w[jpos][kk] +=  buf_ns[kk];

                    }
                } else {
                    /* dealing with lower part */
                    if( jpos == -1 ) {
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
            // new
            for (k=0; k<rnzcount; k++) {
                //~ printf("read %d, %d, %d\n", jrow, k, rnzcount);
                col = rja[k];
                sz = B_DIM(bszc,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       rba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       rba[k], szjrow, zero, buf_ns, dim );
                jpos = iw2[col];

                /* if fill-in element is small then disregard: */
                if( vbnorm2( dim * sz, buf_ns ) < tnorm * drop2 && jpos == -1 ) continue;

                if( jpos == -1 ) {
                    /* this is a fill-in element */
                    jbuf2[len2] = col;
                    iw2[col] = len2;
                    copyBData( dim, sz, w2[len2], buf_ns, 0 );
                    len2++;
                } else {
                    /*-------------------- this is not a fill-in element */
                    for (kk=0;kk<dim*sz; kk++)
                        w2[jpos][kk] +=  buf_ns[kk];
                }
            }
            //end new
            copyBData( dim, szjrow, w[len], buf_fact, 1 );
            jbuf[len] = jrow;
            len++;
        }

        /* update l-matrix */
        lenl = len;
        len = min( lenl, fil0 );
        for( j = 0; j < lenl; j++ ) {
            sz = B_DIM(bsz,jbuf[j]);
            wn[j] = vbnorm2( dim*sz, w[j] );
            iw[j] = j;
        }
        qsplit( wn, iw, &lenl, &len );
        L->nzcount[i] = len;
        if( len > 0 ) {
            L->ja[i] = (int *)Malloc( len * sizeof(int), "vbiluNEW:16" );
            L->ba[i] = (BData *)Malloc( len * sizeof(BData), "vbiluNEW:17" );
        }
        lja = L->ja[i];
        lba = L->ba[i];
        for( j = 0; j < len; j++ ) {
            jpos = iw[j];
            lja[j] = jbuf[jpos];
            sz = B_DIM(bsz,lja[j]);
            lba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:18" );
            //            printf("j value is %d\n", j);//%f %p %s %c

            copyBData( dim, sz, lba[j], w[jpos], 0 );
        }
        for( j = 0; j < lenl; j++ ) iw[j] = -1;

        /* update u-matrix */
        len = min( lenu, fil1 );
        for( j = 1; j < lenu; j++ ) {
            jpos = i+j;
            sz = B_DIM(bsz,jbuf[jpos]);
            wn[j-1] = vbnorm2( dim*sz, w[jpos] );
            iw[j-1] = jpos;
        }
        para = lenu - 1;
        qsplit( wn, iw, &para, &len );
        rnzcount = U->nzcount[i] = len-1;
        if( rnzcount > 0 ) {
            U->ja[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:19" );
            U->ba[i] = (BData *)Malloc( rnzcount*sizeof(BData), "vbiluNEW:20" );
        }
        rja = U->ja[i];
        rba = U->ba[i];
        t = vbnorm2( dim*dim, w[i] );
        for( j = 0; j < rnzcount; j++ ) {
            jpos = iw[j];
            rja[j] = jbuf[jpos];
            sz = B_DIM(bsz,rja[j]);
            rba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:21" );
            copyBData( dim, sz, rba[j], w[jpos], 0 );
            t = max(t,wn[j]);
        }
        for( j = 0; j < lenu-1; j++ ) iw[j] = -1;

        /* save norm in xnrm. Norm = average abs value. divide by len+1
     * instead of len to avoid division by zero. */
        xnrm[i] = t;

        /* store inverse of diagonal element of u */
        D[i] = (BData)Malloc( dim*dim*sizeof(FLOAT), "vbiluNEW:22" );
        copyBData( dim, dim, D[i], w[i], 0 );

#if defined(DBL_CMPLX)
        zgauss( &dim, D[i], &ierr );
#else
        gauss( &dim, D[i], &ierr );
#endif
        if( ierr != 0 ) {
            printf("singular block encountered.\n" );
            for( j = i+1; j < n; j++ ) {
                D[j] = NULL;
                L->ja[j] = NULL;
                L->ba[j] = NULL;
                U->ja[j] = NULL;
                U->ba[j] = NULL;
            }
            return -2;
        }

        for( j = 0; j < lenu; j++ ) {
            iw[ jbuf[i+j] ] = -1;
        }

        /* update l^(-1)f-matrix */
        len = min( len2, fil2 );
        for( j = 0; j < len2; j++ ) {
            sz = B_DIM(bszc,jbuf2[j]);
            wn[j] = vbnorm2( dim*sz, w2[j] );
            iw2[j] = j;
        }
        qsplit( wn, iw2, &len2, &len );
        rnzcount = lflen[i] = len;
        if( rnzcount > 0 ) {
            lfja[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:23" );
            lfma[i] = (BData *)Malloc( rnzcount*sizeof(BData), "vbiluNEW:24" );
        }
        rja = lfja[i];
        rba = lfma[i];
        for( j = 0; j < rnzcount; j++ ) {
            //~ printf("write %d, %d, %d\n", i, j, rnzcount);
            jpos = iw2[j];
            rja[j] = jbuf2[jpos];
            sz = B_DIM(bszc,rja[j]);
            rba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:25" );
            copyBData( dim, sz, rba[j], w2[jpos], 0 );
        }
        for( j = 0; j < len2; j++ ) iw2[j] = -1; //????
        for( j = 0; j < len2; j++ ) iw2[jbuf2[j]] = -1; //????
    }
    //~ for( j = 0; j < n; j++ ) printf("%d ", iw2[j]);
    //~ printf(" done\n");

    bsz = C->bsz;
    bszc = vbmat->E->bszc;
    schur->bsz = (int *)Malloc( (rsize+1)*sizeof(int), "vbiluNEW:26" );
    schur->bsz[rsize] = bsz[rsize];

    /* beginning of second main loop */
    for( i = 0; i < rsize; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        lnzcount = vbmat->E->nzcount[i];
        lja = vbmat->E->ja[i];
        lba = vbmat->E->ba[i];
        rnzcount = C->nzcount[i];
        rja = C->ja[i];
        rba = C->ba[i];

        //~ printf("rnzcount %d, lnzcount %d\n", rnzcount, lnzcount);

        //~ printf("%d %d %d\n", i, rsize, lnzcount);

        tnorm = 0.0;
        for( j = 0; j < rnzcount; j++ ) {
            sz = B_DIM(bsz,rja[j]);
            t = vbnorm2( dim * sz, rba[j] );
            tnorm = max(t,tnorm);
        }
        //~ tnorm = 0.0;

        /* unpack C in arrays w */
        lenl = 0;
        for ( j = 0; j < lnzcount; j++ ) {
            col = lja[j];
            sz = B_DIM(bszc, col);
            jbuf[lenl] = col;
            copyBData( dim, sz, w[lenl], lba[j], 0 );
            iw[col] = lenl;
            lenl++;
        }

        lenu = 0;
        for ( j = 0; j < rnzcount; j++ ) {
            col = rja[j];
            sz = B_DIM(bsz, col);
            jbuf2[lenu] = col;
            copyBData( dim, sz, w2[lenu], rba[j], 0 );
            iw2[col] = lenu;
            lenu++;
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
            szjrow = B_DIM(bszc,jrow);
            if( jpos != j ) {
                col = jbuf[j];
                jbuf[j] = jbuf[jpos];
                jbuf[jpos] = col;
                iw[jrow] = j;
                iw[col] = jpos;
                sz = B_DIM(bszc,col);
                copyBData( dim, sz, buf_ns, w[j], 0 );
                copyBData( dim, szjrow, w[j], w[jpos], 0 );
                copyBData( dim, sz, w[jpos], buf_ns, 0 );
            }
            /* get the multiplier for row to be eliminated (jrow). */
            /* fact = w(jj)*alu(jrow)                              */

            /* NEW */
#if defined(DBL_CMPLX)
            zbxinv( &dim, &szjrow, D[jrow], w[j], buf_fact ); //?????
#else
            bxinv( &dim, &szjrow, D[jrow], w[j], buf_fact ); //?????
#endif             /* zero out element in row by resetting jw(n+jrow) to -1 */
            iw[jrow] = -1;

            if( vbnorm2( dim * szjrow, buf_fact ) * xnrm[jrow] <= tnorm * drop3 )
                continue;

            /* combine current row and row jrow */
            lnzcount = U->nzcount[jrow];
            lja = U->ja[jrow];
            lba = U->ba[jrow];
            rnzcount = lflen[jrow];
            rja = lfja[jrow];
            rba = lfma[jrow];
            for( k = 0; k < lnzcount; k++ ) {
                col = lja[k];
                sz = B_DIM(bszc,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       lba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       lba[k], szjrow, zero, buf_ns, dim );
                jpos = iw[col];

                /* if fill-in element is small then disregard: */
                if( vbnorm2( dim * sz, buf_ns ) < tnorm * drop4 && jpos == -1 ) continue;

                if( jpos == -1 ) {
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
            // new
            for (k=0; k<rnzcount; k++) {
                col = rja[k];
                sz = B_DIM(bsz,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       rba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       rba[k], szjrow, zero, buf_ns, dim );
                jpos = iw2[col];

                /* if fill-in element is small then disregard: */
                if( vbnorm2( dim * sz, buf_ns ) < tnorm * drop4 && jpos == -1 ) continue;

                if( jpos == -1 ) {
                    /* this is a fill-in element */
                    jbuf2[lenu] = col;
                    iw2[col] = lenu;
                    copyBData( dim, sz, w2[lenu], buf_ns, 0 );
                    lenu++;
                } else {
                    /*-------------------- this is not a fill-in element */
                    for (kk=0;kk<dim*sz; kk++)
                        w2[jpos][kk] +=  buf_ns[kk];
                }
            }
            //end new
            copyBData( dim, szjrow, w[len], buf_fact, 1 );
            jbuf[len] = jrow;
            len++;
        }

        if( lenu == 0 ) {
            jbuf2[0] = i;
            iw2[0] = 0;
            for (j = 0; j < dim; j++)
                for (k = 0; k < dim; k++)
                    w2[0][j*dim+k] = (j == k) ? 1.0 : 0.0;
            //~ copyBData( dim, sz, w2[0], lba[j], 0 );
            lenu = 1;
        }

        /* update l-matrix */
        for( j = 0; j < lenl; j++ ) iw[j] = -1;

        /* update schur-matrix */
        len = min( lenu, fil4 );
        for( j = 0; j < lenu; j++ ) {
            sz = B_DIM(bsz,jbuf2[j]);
            wn[j] = vbnorm2( dim*sz, w2[j] );
            iw2[j] = j;
        }
        qsplit( wn, iw2, &lenu, &len );
        rnzcount = schur->nzcount[i] = len;
        //~ printf("%d\n", rnzcount);
        if( rnzcount > 0 ) {
            schur->ja[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:27" );
            schur->ba[i] = (BData *)Malloc( rnzcount*sizeof(BData), "vbiluNEW:28" );
        }
        rja = schur->ja[i];
        rba = schur->ba[i];
        schur->bsz[i] = bsz[i];
        //~ t = wn[0];
        for( j = 0; j < rnzcount; j++ ) {
            jpos = iw2[j];
            rja[j] = jbuf2[jpos];
            sz = B_DIM(bsz,rja[j]);
            rba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:29" );
            copyBData( dim, sz, rba[j], w2[jpos], 0 );
            //~ t = max(t,wn[j]);
        }
        for( j = 0; j < lenu; j++ ) iw2[j] = -1;
        for( j = 0; j < lenu; j++ ) iw2[jbuf2[j]] = -1;

        /* save norm in xnrm. Norm = average abs value. divide by len+1
     * instead of len to avoid division by zero. */
        //~ xnrm[i] = t; //?????


        //~ D[i] = (BData)Malloc( dim*dim*sizeof(FLOAT), "vbiluNEW" );
        //~ copyBData( dim, dim, D[i], w[i], 0 );
        //~ gauss( &dim, D[i], &ierr );

        //~ printf("%d, %d, %d, %d\n", i, lenu, lenl, rsize);
    }

    vbmat->lu->DiagOpt = 1;
    free( jbuf );
    free( jbuf2 );
    free( buf_fact );
    free( buf_ns );
    free( xnrm );
    free( wn );
    free( iw );
    free( iw2 );
    for (i=0; i<lsize; i++) {
        if (lflen[i] > 0) {
            for (j=0; j<lflen[i]; j++)
                free(lfma[i][j]);
            free(lfma[i]);
            free(lfja[i]);
        }
    }
    for (i=0; i<n; i++) {
        free(w[i]);
        free(w2[i]);
    }
    free( w );
    free( w2 );
    free(lfma);
    free(lfja);
    free(lflen);

    return 0;
}

int vbilukNEW(vbp4ptr vbmat, vbsptr B, vbsptr C, 
              int lfil, vbsptr schur)
{
    /*----------------------------------------------------------------------------
 * Block ILUT (BILUT) preconditioner
 * Block incomplete LU factorization with dual truncation mechanism
 * NOTE : no pivoting implemented as yet in GE for diagonal blocks
 *----------------------------------------------------------------------------
 * Parameters
 *----------------------------------------------------------------------------
 * on entry:
 * =========
 * vbmat    = block matrix stored in VBSpaFmt format -- see globheads.h for
 *            details on format, the block sizes might be different
 * lu       = pointer to a VBILUSpar struct -- see globheads.h for details
 *            on format
 * lfil     = integer. The fill-in parameter. Each row of L and
 *            each row of U will have a maximum of lfil elements.
 *            WARNING: THE MEANING OF LFIL HAS CHANGED WITH RESPECT TO
 *            EARLIER VERSIONS.
 *            lfil must be .ge. 0.
 * tol      = real*8. Sets the threshold for dropping small terms in the
 *            factorization. See below for details on dropping strategy.
 * w        = working array
 * fp       = file pointer for error log ( might be stdout )
 *
 * on return:
 * ==========
 * ierr     = return value.
 *            ierr  = 0   --> successful return.
 *            ierr  = -1  --> Illegal value for lfil
 *            ierr  = -2  --> singular block or zero row encountered
 * lu->n    = dimension of the block matrix
 *   ->bsz  = the row/col of the first element of each diagonal block
 *            the size of the i-th row block should be bsz[i+1] - bsz[i]
 *   ->L    = L part -- stored in VBSpaFmt format
 *   ->D    = Diagonals
 *   ->U    = U part -- stored in VBSpaFmt format
 *----------------------------------------------------------------------------
 * Notes:
 * ======
 * All the diagonal blocks of the input block matrix must not be singular
 *----------------------------------------------------------------------------
 * Dual drop-off strategy works as follows.
 *
 * 1) Theresholding in L and U as set by tol. Any element whose size
 *    is less than some tolerance (relative to the norm of current
 *    row in u) is dropped.
 *
 * 2) Keeping only the largest lfil elements in the i-th row of L
 *    and the largest lfil elements in the i-th row of U.
 *
 * Flexibility: one can use tol=0 to get a strategy based on keeping the
 * largest elements in each row of L and U. Taking tol .ne. 0 but lfil=n
 * will give the usual threshold strategy (however, fill-in is then
 * impredictible).
 *--------------------------------------------------------------------------*/
    int n = schur->n, *bsz = vbmat->bsz, *bszc = vbmat->F->bszc, ierr;
    FLOAT one = 1.0, zero = 0.0;
    int dim, szjrow, sz, len, lenu, lenl, len2, col, jpos, jrow, upos;
    int max_blk_sz = MAX_BLOCK_SIZE*MAX_BLOCK_SIZE*sizeof(double);
    int lnzcount, rnzcount, *jbuf, *iw, i, j, k, kk, *lja, *rja;
    int *jbuf2, *iw2, lsize, rsize, **lfja, *lflen;
    int *levls, *levls2, **levs, **levs2, *lev, jlev;
    vbsptr L, U;
    BData *lba, *rba, *D, buf_fact, buf_ns;
    BData *w, *w2, **lfma;

    lsize = vbmat->nB;
    rsize = C->n;
    printf("lsize %d, rsize %d, n %d\n", lsize, rsize, n);
    n = lsize > rsize ? lsize : rsize;

    if(lfil < 0 || n <= 0) {
        printf( "vbiluNEW: Illegal value for lfil.\n" );
        return -1;
    }
    buf_fact = (BData)Malloc( max_blk_sz, "vbiluNEW:1" );
    buf_ns = (BData)Malloc( max_blk_sz, "vbiluNEW:2" );
    levls = (int *)Malloc( n * sizeof(int), "vbiluNEW:3" );
    levs = (int **)Malloc( n * sizeof(int *), "vbiluNEW:4" );
    levls2 = (int *)Malloc( n * sizeof(int), "vbiluNEW:5" );
    levs2 = (int **)Malloc( n * sizeof(int *), "vbiluNEW:6" );

    jbuf = (int *)Malloc( n * sizeof(int), "vbiluNEW:5" );
    w = (BData *)Malloc( n * sizeof(BData),"vbiluNEW:6");
    for( i = 0; i < n; i++ )
        w[i] = (FLOAT *)Malloc( max_blk_sz, "vbiluNEW:7" );
    iw = (int *)Malloc( n * sizeof(int), "vbiluNEW:8" );

    L = vbmat->lu->L;
    U = vbmat->lu->U;
    D = vbmat->lu->D;

    jbuf2 = (int *)Malloc( n * sizeof(int), "vbiluNEW:9" );
    w2 = (BData *)Malloc( n * sizeof(BData),"vbiluNEW:10");
    for( i = 0; i < n; i++ )
        w2[i] = (FLOAT *)Malloc( max_blk_sz, "vbiluNEW:11" );
    iw2 = (int *)Malloc( n * sizeof(int), "vbiluNEW:12" );

    lfma = (BData **) Malloc(lsize*sizeof(BData *), "vbiluNEW:13" );
    lfja = (int **) Malloc(lsize*sizeof(int *), "vbiluNEW:14" );
    lflen = (int *) Malloc(lsize*sizeof(int), "vbiluNEW:15" );

    /* set indicator array jw to -1 */
    for( i = 0; i < n; i++ )
    {
        iw[i] = -1;
        iw2[i] = -1;
        levls[i] = 0;
        levls2[i] = 0;
    }

    /* beginning of main loop */
    for( i = 0; i < lsize; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        lnzcount = B->nzcount[i];
        lja = B->ja[i];
        lba = B->ba[i];
        rnzcount = vbmat->F->nzcount[i];
        rja = vbmat->F->ja[i];
        rba = vbmat->F->ba[i];

        /* unpack L-part and U-part of row of A in arrays w */
        lenu = 1;
        lenl = 0;
        jbuf[i] = i;
        zrmC( dim, dim, w[i] );
        iw[i] = i;

        for( j = 0; j < lnzcount; j++ ) {
            col = lja[j];
            sz = B_DIM(bsz, col);
            if( col < i ) {
                jbuf[lenl] = col;
                levls[lenl] = 0;
                iw[col] = lenl;
                copyBData( dim, sz, w[lenl], lba[j], 0 );
                lenl++;
            } else if( col == i ) {
                levls[i] = 0;
                copyBData( dim, dim, w[i], lba[j], 0 );
            } else {
                jpos = i + lenu;
                jbuf[jpos] = col;
                levls[jpos] = 0;
                iw[col] = jpos;
                copyBData( dim, sz, w[jpos], lba[j], 0 );
                lenu++;
            }
        }

        len2 = 0;
        for ( j = 0; j < rnzcount; j++ ) {
            col = rja[j];
            sz = B_DIM(bszc, col);
            jbuf2[len2] = col;
            levls2[len2] = 0;
            copyBData( dim, sz, w2[len2], rba[j], 0 );
            iw2[col] = len2;
            len2++;
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
#endif                /* zero out element in row by resetting jw(n+jrow) to -1 */
            iw[jrow] = -1;

            lev = levs[jrow];
            jlev = levls[j];

            if (jlev > lfil) continue;

            /* combine current row and row jrow */
            lnzcount = U->nzcount[jrow];
            lja = U->ja[jrow];
            lba = U->ba[jrow];

            rnzcount = lflen[jrow];
            rja = lfja[jrow];
            rba = lfma[jrow];

            for( k = 0; k < lnzcount; k++ ) {
                col = lja[k];
                sz = B_DIM(bsz,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       lba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       lba[k], szjrow, zero, buf_ns, dim );
                jpos = iw[col];

                if( col >= i ) {
                    /* dealing with upper part */
                    //          if( jpos == -1 ) {
                    if( jpos == -1) {
                        /* this is a fill-in element */
                        upos = i + lenu;
                        jbuf[upos] = col;
                        levls[upos] = jlev+lev[k]+1;
                        iw[col] = upos;
                        copyBData( dim, sz, w[upos], buf_ns, 0 );
                        lenu++;
                    } else {
                        /*-------------------- this is not a fill-in element */
                        for (kk=0;kk<dim*sz; kk++)
                            w[jpos][kk] +=  buf_ns[kk];
                        levls[jpos] = min(levls[jpos], jlev+lev[k]+1);

                    }
                } else {
                    /* dealing with lower part */
                    if( jpos == -1 ) {
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

            // new
            lev = levs2[jrow];
            jlev = levls2[j];

            if (jlev > lfil) continue;

            for (k=0; k<rnzcount; k++) {
                //~ printf("read %d, %d, %d\n", jrow, k, rnzcount);
                col = rja[k];
                sz = B_DIM(bszc,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       rba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       rba[k], szjrow, zero, buf_ns, dim );
                jpos = iw2[col];

                if( jpos == -1 ) {
                    /* this is a fill-in element */
                    jbuf2[len2] = col;
                    levls2[len2] = jlev+lev[k]+1;
                    iw2[col] = len2;
                    copyBData( dim, sz, w2[len2], buf_ns, 0 );
                    len2++;
                } else {
                    /*-------------------- this is not a fill-in element */
                    for (kk=0;kk<dim*sz; kk++)
                        w2[jpos][kk] +=  buf_ns[kk];
                    levls2[jpos] = min(levls2[jpos], jlev+lev[k]+1);
                }
            }
            //end new
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
        L->nzcount[i] = len;
        if( len > 0 ) {
            L->ja[i] = (int *)Malloc( len * sizeof(int), "vbiluNEW:16" );
            L->ba[i] = (BData *)Malloc( len * sizeof(BData), "vbiluNEW:17" );
        }
        lja = L->ja[i];
        lba = L->ba[i];
        for( j = 0; j < len; j++ ) {
            jpos = iw[j];
            lja[j] = jbuf[jpos];
            sz = B_DIM(bsz,lja[j]);
            lba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:18" );
            //            printf("j value is %d\n", j);//%f %p %s %c

            copyBData( dim, sz, lba[j], w[jpos], 0 );
        }
        for( j = 0; j < lenl; j++ ) iw[j] = -1;

        /* update u-matrix */
        len = 0;
        for( j = 1; j < lenu; j++ ) {
            jpos = i+j;
            if (levls[jpos] <= lfil)
            {
                iw[len] = jpos;
                len++;
            }
        }
        rnzcount = U->nzcount[i] = len;
        if( rnzcount > 0 ) {
            U->ja[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:19" );
            U->ba[i] = (BData *)Malloc( rnzcount*sizeof(BData), "vbiluNEW:20" );
            levs[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:21" );
        }
        rja = U->ja[i];
        rba = U->ba[i];
        lev = levs[i];

        for( j = 0; j < rnzcount; j++ ) {
            jpos = iw[j];
            rja[j] = jbuf[jpos];
            lev[j] = levls[jpos];
            sz = B_DIM(bsz,rja[j]);
            rba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:21" );
            copyBData( dim, sz, rba[j], w[jpos], 0 );
        }
        for( j = 0; j < lenu; j++ ) iw[j] = -1;

        /* store inverse of diagonal element of u */
        D[i] = (BData)Malloc( dim*dim*sizeof(FLOAT), "vbiluNEW:22" );
        copyBData( dim, dim, D[i], w[i], 0 );

#if defined(DBL_CMPLX)
        zgauss( &dim, D[i], &ierr );
#else
        gauss( &dim, D[i], &ierr );
#endif
        if( ierr != 0 ) {
            printf("singular block encountered.\n" );
            for( j = i+1; j < n; j++ ) {
                D[j] = NULL;
                L->ja[j] = NULL;
                L->ba[j] = NULL;
                U->ja[j] = NULL;
                U->ba[j] = NULL;
            }
            return -2;
        }

        for( j = 0; j < lenu; j++ ) {
            iw[ jbuf[i+j] ] = -1;
        }

        /* update l^(-1)f-matrix */
        len = 0;
        for( j = 0; j < len2; j++ ) {
            if (levls2[j] <= lfil) {
                iw2[len] = j;
                len++;
            }
        }
        rnzcount = lflen[i] = len;
        if( rnzcount > 0 ) {
            lfja[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:23" );
            lfma[i] = (BData *)Malloc( rnzcount*sizeof(BData), "vbiluNEW:24" );
            levs2[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:25" );
        }
        rja = lfja[i];
        rba = lfma[i];
        lev = levs2[i];
        for( j = 0; j < rnzcount; j++ ) {
            //~ printf("write %d, %d, %d\n", i, j, rnzcount);
            jpos = iw2[j];
            rja[j] = jbuf2[jpos];
            lev[j] = levls2[jpos];
            sz = B_DIM(bszc,rja[j]);
            rba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:25" );
            copyBData( dim, sz, rba[j], w2[jpos], 0 );
        }
        for( j = 0; j < len2; j++ ) iw2[j] = -1; //????
        for( j = 0; j < len2; j++ ) iw2[jbuf2[j]] = -1; //????
    }
    //~ for( j = 0; j < n; j++ ) printf("%d ", iw2[j]);
    //~ printf(" done\n");

    bsz = C->bsz;
    bszc = vbmat->E->bszc;
    schur->bsz = (int *)Malloc( (rsize+1)*sizeof(int), "vbiluNEW:26" );
    schur->bsz[rsize] = bsz[rsize];

    /* beginning of second main loop */
    for( i = 0; i < rsize; i++ ) {
        dim = B_DIM(bsz,i);  /* number of rows of blocks in i-th row */
        lnzcount = vbmat->E->nzcount[i];
        lja = vbmat->E->ja[i];
        lba = vbmat->E->ba[i];
        rnzcount = C->nzcount[i];
        rja = C->ja[i];
        rba = C->ba[i];

        /* unpack C in arrays w */
        lenl = 0;
        for ( j = 0; j < lnzcount; j++ ) {
            col = lja[j];
            sz = B_DIM(bszc, col);
            jbuf[lenl] = col;
            levls[lenl] = 0;
            copyBData( dim, sz, w[lenl], lba[j], 0 );
            iw[col] = lenl;
            lenl++;
        }

        lenu = 0;
        for ( j = 0; j < rnzcount; j++ ) {
            col = rja[j];
            sz = B_DIM(bsz, col);
            jbuf2[lenu] = col;
            levls2[lenu] = 0;
            copyBData( dim, sz, w2[lenu], rba[j], 0 );
            iw2[col] = lenu;
            lenu++;
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
            szjrow = B_DIM(bszc,jrow);
            if( jpos != j ) {
                col = jbuf[j];
                jbuf[j] = jbuf[jpos];
                jbuf[jpos] = col;
                iw[jrow] = j;
                iw[col] = jpos;
                sz = B_DIM(bszc,col);
                col = levls[j];
                levls[j] = levls[jpos];
                levls[jpos] = col;
                copyBData( dim, sz, buf_ns, w[j], 0 );
                copyBData( dim, szjrow, w[j], w[jpos], 0 );
                copyBData( dim, sz, w[jpos], buf_ns, 0 );
            }
            /* get the multiplier for row to be eliminated (jrow). */
            /* fact = w(jj)*alu(jrow)                              */

            /* NEW */
#if defined(DBL_CMPLX)
            zbxinv( &dim, &szjrow, D[jrow], w[j], buf_fact ); //?????
#else
            bxinv( &dim, &szjrow, D[jrow], w[j], buf_fact ); //?????
#endif             /* zero out element in row by resetting jw(n+jrow) to -1 */
            iw[jrow] = -1;

            lev = levs[jrow];
            jlev = levls[j];

            if (jlev > lfil) continue;

            /* combine current row and row jrow */
            lnzcount = U->nzcount[jrow];
            lja = U->ja[jrow];
            lba = U->ba[jrow];
            rnzcount = lflen[jrow];
            rja = lfja[jrow];
            rba = lfma[jrow];
            for( k = 0; k < lnzcount; k++ ) {
                col = lja[k];
                sz = B_DIM(bszc,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       lba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       lba[k], szjrow, zero, buf_ns, dim );
                jpos = iw[col];

                if( jpos == -1 ) {
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
            // new
            lev = levs2[jrow];
            jlev = levls2[j];

            if (jlev > lfil) continue;

            for (k=0; k<rnzcount; k++) {
                col = rja[k];
                sz = B_DIM(bsz,col);
                //                dgemm ("n","n",&dim, &sz, &szjrow, &one, buf_fact, &dim,
                //                       rba[k], &szjrow, &zero, buf_ns, &dim );
                GGEMM ("n", "n", dim, sz, szjrow, one, buf_fact, dim,
                       rba[k], szjrow, zero, buf_ns, dim );
                jpos = iw2[col];

                if( jpos == -1 ) {
                    /* this is a fill-in element */
                    jbuf2[lenu] = col;
                    levls2[lenu] = jlev+lev[k]+1;
                    iw2[col] = lenu;
                    copyBData( dim, sz, w2[lenu], buf_ns, 0 );
                    lenu++;
                } else {
                    /*-------------------- this is not a fill-in element */
                    for (kk=0;kk<dim*sz; kk++)
                        w2[jpos][kk] +=  buf_ns[kk];
                    levls2[jpos] = min(levls2[jpos], jlev+lev[k]+1);
                }
            }
            //end new
            copyBData( dim, szjrow, w[j], buf_fact, 1 );
            jbuf[j] = jrow;
        }

        if( lenu == 0 ) {
            jbuf2[0] = i;
            iw2[0] = 0;
            for (j = 0; j < dim; j++)
                for (k = 0; k < dim; k++)
                    w2[0][j*dim+k] = (j == k) ? 1.0 : 0.0;
            lenu = 1;
        }

        /* update l-matrix */
        for( j = 0; j < lenl; j++ ) iw[j] = -1;

        /* update schur-matrix */
        len = 0;
        for( j = 0; j < lenu; j++ ) {
            if (levls2[j] <= lfil) {
                iw2[len] = j;
                len++;
            }
        }
        rnzcount = schur->nzcount[i] = len;
        //~ printf("%d\n", rnzcount);
        if( rnzcount > 0 ) {
            schur->ja[i] = (int *)Malloc( rnzcount*sizeof(int), "vbiluNEW:27" );
            schur->ba[i] = (BData *)Malloc( rnzcount*sizeof(BData), "vbiluNEW:28" );
        }
        rja = schur->ja[i];
        rba = schur->ba[i];
        schur->bsz[i] = bsz[i];
        //~ t = wn[0];
        for( j = 0; j < rnzcount; j++ ) {
            jpos = iw2[j];
            rja[j] = jbuf2[jpos];
            sz = B_DIM(bsz,rja[j]);
            rba[j] = (BData)Malloc( dim*sz*sizeof(FLOAT), "vbiluNEW:29" );
            copyBData( dim, sz, rba[j], w2[jpos], 0 );
        }
        for( j = 0; j < lenu; j++ ) iw2[j] = -1;
        for( j = 0; j < lenu; j++ ) iw2[jbuf2[j]] = -1;
    }

    vbmat->lu->DiagOpt = 1;
    free( jbuf );
    free( jbuf2 );
    free( buf_fact );
    free( buf_ns );
    free( iw );
    free( iw2 );
    for (i=0; i<lsize; i++) {
        if (lflen[i] > 0) {
            for (j=0; j<lflen[i]; j++)
                free(lfma[i][j]);
            free(lfma[i]);
            free(lfja[i]);
            free(levs2[i]);
        }
        if (U->nzcount[i] > 0)
            free(levs[i]);
    }
    for (i=0; i<n; i++) {
        free(w[i]);
        free(w2[i]);
    }
    free( w );
    free( w2 );
    free(lfma);
    free(lfja);
    free(lflen);
    free(levls);
    free(levs);
    free(levls2);
    free(levs2);

    return 0;
}
