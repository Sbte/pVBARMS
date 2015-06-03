#ifndef __ITSOL_INCLUDED_PROTOS_H__
#define __ITSOL_INCLUDED_PROTOS_H__

#include "globheads.h"

#if defined(FORTRAN_CAPS)
#define qsplit  QSPLIT
#define readmtc READMTC
#define csrcsc  CSRCSC
#define roscal  ROSCAL
#define coscal  COSCAL
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define qsplit  qsplit__
#define readmtc readmtc__
#define csrcsc  csrcsc__
#define roscal  roscal__
#define coscal  coscal__
#elif defined(FORTRAN_UNDERSCORE)
#define qsplit  qsplit_
#define readmtc readmtc_
#define csrcsc  csrcsc_
#define roscal  roscal_
#define coscal  coscal_
#else
#include <essl.h>
#endif

#ifndef min
#define min(a,b) (((a)>(b))?(b):(a))
#endif
#ifndef max
#define max(a,b) (((a)>(b))?(a):(b))
#endif

/* sets */
extern void errexit(char *f_str, ...);
extern void *Malloc(int nbytes, const char *msg); 
void *Calloc( int n, int size, const char *msg );
extern int setupP4 (p4ptr amat, int Bn, int Cn,  csptr F,  csptr E);
extern int cleanP4(p4ptr amat);
extern int setupILUT(ilutptr amat, int len);
extern int cleanILUT(ilutptr amat, int indic);
extern int setupCS(csptr amat, int len, int job); 
extern int cleanCS(csptr amat);
extern int nnzCS(csptr amat); 
//extern int outputcsmat ( csptr mat, char *filename, int onebase);//new
extern int nnz_arms (arms PreSt);
extern int cs_nnz (csptr A) ;
extern int cscpy(csptr amat, csptr bmat);
extern int cleanILU( iluptr lu );
extern int mallocRow( iluptr lu, int nrow );
extern int CSRcs(int n, FLOAT *a, int *ja, int *ia, csptr mat, int rsa);
extern int lev4_nnz(p4ptr levmat, int *lev); 
extern int setupILU( iluptr lu, int n );
extern void setup_arms (arms Levmat);
extern int cleanARMS(arms ArmsPre);
extern int csSplit4(csptr amat, int bsize, int csize, csptr B, csptr F,
		    csptr E, csptr C);
//----------------------------------------------------------------dividing line-----------------------------------------------------------------------
extern int setupVBMat( vbsptr vbmat, int n, int *nB );
extern int setupVBMat1( vbsptr vbmat, int n,int nc, int *nBR, int *nBC );
extern int col2vbcol( int col, vbsptr vbmat );
extern int csrvbsrC( int job, int nBlk, int *nB, csptr csmat, vbsptr vbmat );
extern int csrvbsrC_new( int job, int nBlk, int *nB, csptr csmat, vbsptr vbmat );
extern int nnzVBMat( vbsptr vbmat );
extern int nnzVBMat1( vbsptr vbmat );
extern int memVBMat( vbsptr vbmat );
extern int memVBMat1( vbsptr vbmat );
extern int mallocVBRow( vbiluptr lu, int nrow );
extern void zrmC( int m, int n, BData data );
extern void copyBData( int m, int n, BData dst, BData src, int isig );
extern int permuten(int *iwork, int *iworkn, int *bsz, int n);
//----------------------------------------------------------------dividing line-----------------------------------------------------------------------
/* MatOps */
extern void matvec(csptr mata, FLOAT *x, FLOAT *y); 
extern void matvecz(csptr mata, FLOAT *x, FLOAT *y, FLOAT *z);
extern void luinv(int n, FLOAT *a, FLOAT *x, FLOAT *y); 
extern int lusolC( FLOAT *y, FLOAT *x, iluptr lu ); 
extern int rpermC(csptr mat, int *perm); 
extern int cpermC(csptr mat, int *perm) ; 
extern int dpermC(csptr mat, int *perm) ; 
extern int CSparTran(csptr amat, csptr bmat, CompressType *compress);
extern void invsp(int start, ilutptr ilusch, FLOAT *y, FLOAT *x);
extern int armsol2(FLOAT *x,  arms Prec);
extern int ascend (p4ptr levmat, FLOAT *x, FLOAT *wk);
extern int descend(p4ptr levmat, FLOAT *x, FLOAT *wk);
extern p4ptr Lvsol2(FLOAT *x, int nlev, p4ptr levmat, ilutptr ilusch,
		    int flag); 
extern int   Uvsol2(FLOAT *x, int nlev, int n, p4ptr levmat, ilutptr
		    ilusch); 
extern void  SchLsol(ilutptr ilusch, FLOAT *y) ;
extern void  SchUsol(ilutptr ilusch, FLOAT *y) ;
extern void Lsolp(int start, csptr mata, FLOAT *b, FLOAT *y);
extern void Usolp(int start, csptr mata, FLOAT *y, FLOAT *x);
extern int invGauss(int nn, FLOAT *A); 
extern int invSVD(int nn, FLOAT *A) ;
extern void arms_Usol(csptr mata, FLOAT *b, FLOAT *x);
extern void arms_Lsol(csptr mata, FLOAT *b, FLOAT *x);
extern int condestArms(arms armspre, FLOAT *y, FILE *fp );

/* misc.c */
extern int SparTran(csptr amat, csptr bmat, int job, int flag); 
extern int coscalC(csptr mata, double *diag, int nrm);
extern void dscale(int n, double *dd, FLOAT *x, FLOAT * y);
extern void hilosort(csptr mat, int abval, int hilo);
extern void printmat(FILE *ft, csptr A, int i0, int i1);
extern void qqsort(int *ja, FLOAT *ma, int left, int right);
extern void qsort2C(int *ja, FLOAT *ma, int left, int right, int
		    abval); 
extern void qsort3i(int *wa, int *cor1, int *cor2, int left, int
		    right); 
extern void qsortC(int *ja, FLOAT *ma, int left, int right, int
		   abval); 
extern void qsortR2I(double *wa, int *cor1, int *cor2, int left, int
		     right); 
extern int qsplitC(FLOAT *a, int *ind, int n, int ncut);
extern int roscalC(csptr mata, double *diag, int nrm);
extern void swapj(int v[], int i, int j);
extern void swapm(FLOAT v[], int i, int j);
extern void dswapm(double v[], int i, int j);



/* PQ.c */
extern int PQperm(csptr mat, int *Pord, int *Qord, int
		  *nnod, double tol, int nbnd);
extern int preSel(csptr mat, int *icor, int *jcor, int job, double
		  tol, int *count, int nbnd);
extern int weightsC(csptr mat, double *w);
extern int add2com(int *nback, int nod, int *iord, int *riord);
extern int add2is(int *last, int nod, int *iord, int *riord);
extern int indsetC(csptr mat, int bsize, int *iord, int *nnod, double
		   tol,int nbnd); 

/* setblks.c */
extern int KeyComp( const void *vfst, const void *vsnd );
extern int init_blocks( csptr csmat, int *pnBlock, int **pnB, int
			**pperm, double eps);//, double *t_hash, double*t_angle );
/*cliques.c*/
extern int init_blocks_density(csptr A, int *block_number, int **block_sizes, int **pperm, double eps);

/* systimer.c */
extern double sys_timer();

/* barmsmisc.c */
extern int vbsrc2csr(vbsptr vbmat, csptr csmat);
extern int vbrpermC(vbsptr vbmat, int *perm);
extern int vbcpermC(vbsptr vbmat, int *perm); 
extern int vbSplit4copy(vbsptr vbmat, int bsize, int csize, vbsptr B, vbsptr F, vbsptr E, vbsptr C);
extern int setupVBP4 (vbp4ptr vbmat, int Bn, int Cn,  vbsptr F,  vbsptr E , int *bsz);
extern int setupVBP4_vbarmsold (vbp4ptr vbmat, int Bn, int Cn,  vbsptr F,  vbsptr E , int *bsz);
extern int vbilutD( vbsptr schur, double *droptol, int *lfil, vbilutptr vbmat);
extern int setupVBILUT(vbilutptr vbmat, int len, int *bsz);
extern int permuten(int *iwork, int *iworkn, int *bsz, int n);
extern int vbascend (vbp4ptr levmat, FLOAT *x, FLOAT *wk);
extern void VBSchLUsol(vbilutptr ilusch, FLOAT *y);
extern void pSchLUsol(vbilutptr ilusch, FLOAT *y);
extern int barmsol2(FLOAT *x,  vbarms Prec);
extern int cleanBARMS(vbarms ArmsPre);
extern int cleanVBMat(vbsptr mata);
extern int nnz_vbarms (vbarms PreSt,  FILE *ft);
extern void vbdscale(int nn, double *dd, FLOAT *x, FLOAT * y);
extern int vbfroscalC(vbsptr mata, double *diag, int nrm);
extern int vbfcoscalC(vbsptr mata, double *diag, int nrm);
extern int indsetC2(csptr mat, int bsize, int *iord, int *nnod, double tol, int nbnd, int *nBB, int *nset);
extern vbp4ptr BLvsol2(FLOAT *x, int nlev, vbp4ptr levmat, vbilutptr ilusch, int flag);
extern void vbLsolp(int start, vbiluptr lu, FLOAT *y, FLOAT *x);
extern void vbUsolp(int start, vbiluptr lu, FLOAT *y, FLOAT *x);
extern void vbinvsp(int start, vbiluptr lu, FLOAT *y, FLOAT *x);
extern int vblusolC( FLOAT *y, FLOAT *x, vbiluptr lu);
extern int cleanVBILU( vbiluptr lu );
extern int setupVBILU( vbiluptr lu, int n, int *bsz );
extern int nnz_vbilu(vbiluptr lu );
extern int vbiluNEW(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol,int *lfil, vbsptr schur);
extern int vbilukNEW(vbp4ptr vbmat, vbsptr B, vbsptr C, int lfil, vbsptr schur);
extern int computschurpartition(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol, int *lfil, vbsptr schur, int *nBB, int nset);
extern int vbilutC( vbsptr vbmat, vbiluptr lu, int lfil, double tol, BData *w, FILE *fp );
extern int vbilukC( int lofM, vbsptr vbmat, vbiluptr lu, FILE *fp );
extern int computLUFsparse(vbsptr F, vbiluptr lu, vbsptr FF);
extern int cleanVBMat1( vbsptr vbmat );
//~ extern int computschur(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol, int *lfil, vbsptr schur);
//~ extern int vbplussvb(int job, vbsptr A, double s, vbsptr B, vbsptr C);
//~ extern int vbsrc2csr1(vbsptr vbmat, csptr csmat);
//~ extern int outputvbmat( vbsptr vbmat, char *filename, int onebase);
//~ extern int outputvbmat1( vbsptr vbmat, char *filename, int onebase);
//~ extern int vbSplit4(vbsptr vbmat, int bsize, int csize, vbsptr B, vbsptr F, vbsptr E, vbsptr C);;
//~ extern int outputBData( BData *mat, char *filename, int onebase,int n,int nc);
//~ extern int vb2BData(vbsptr vbmat, BData *mat);
//~ extern int BData2vb(BData *mat, vbsptr vbmat, int n, int nc, int *bsz, int *bszc);
//~ extern int outputvbD( vbiluptr lu, char *filename, int onebase);
//~ extern int vbplussvb1(int job, vbsptr A, double s, vbsptr B, vbsptr C);//failed
//~ extern int computLUF(vbsptr F, vbiluptr lu, vbsptr FF);
//~ extern int vbmulvb(int job, vbsptr A, vbsptr B, vbsptr C);
//~ extern void vbtest_ssop();
//~ extern int vbscpy(vbsptr amat, vbsptr bmat);
//~ extern void permuterhs(double *rhs, double *rhs1, int nn, int *perm);
//~ extern void vbmatvecz(vbsptr vbmat, double *x, double *y, double *z );
//~ extern int vbdescend(vbp4ptr levmat, double *x, double *wk);
//~ extern void vbmatvec1(vbsptr vbmat, double *x, double *y );
//~ extern void permuterhs(double *rhs, double *rhs1, int nn, int *perm);
//~ extern void vbmatvecz(vbsptr vbmat, double *x, double *y, double *z );
//~ extern int vbdescend(vbp4ptr levmat, double *x, double *wk);
//~ extern void vbmatvec1(vbsptr vbmat, double *x, double *y );
//~ extern vbp4ptr BLvsol2(double *x, int nlev, vbp4ptr levmat, vbilutptr ilusch);
//~ extern int BUvsol2(double *x, int nlev, int n, vbp4ptr levmat, vbilutptr ilusch);
//~ extern int condestBArms(vbarms armspre, double *y, FILE *fp );
//~ extern int cleanVBILUT(vbilutptr amat, int indic);
//~ extern int cleanVBP4(vbp4ptr amat);
//~ extern int vblev4_nnz(vbp4ptr levmat, int *lev, FILE *ft);
//~ extern int vbplussvb1(int job, vbsptr A, double s, vbsptr B, vbsptr C, double tol);
//~ extern int vbroscalC(vbsptr mata, double *diag, int nrm);
//~ extern int vbcoscalC1(vbsptr mata, double *diag, int nrm);
//~ extern void generatefull(double *d1, double *d11, int n, int *bsz);
//~ extern int vbcoscalC(vbsptr mata, double *diag, int nrm);
//~ extern void vproduct(int n, double *dst, double *src, int isig);
//~ extern int setupVBvector( vbbsptr vbbmat, int nset, int *nBB);
//~ extern int setupVBilutvector( vbbiluptr vbbmat, int nset, int *nBB);
//~ extern int Bsplit(vbsptr B, vbbsptr Bb, vbsptr F, vbbsptr Ff, int *nBB, int nset);
//~ extern int computLUFinpartition( vbbsptr Bb, vbbsptr Ff, vbbsptr FF, vbbiluptr luu,io_t io, FILE *flog);
//~ extern int Ffassemble( vbsptr FF, vbbsptr Ff);
//~ extern int BLUassemble( vbiluptr lu, vbbiluptr luu);
//~ extern int cleanVBBMat( vbbsptr vbbmat );
//~ extern int cleanVBBILU( vbbiluptr luu );
//~ extern int armsmb2(csptr Amat, int *ipar, double *droptol, int *lfil, double tolind, arms PreMat, FILE *ft);
//~ extern void setup_mvbarms (mvbarms Levmat);
//~ extern int MBarms2(csptr Amat, int *ipar, double *droptol, int *lfil, double tolind, mvbarms PreMatm, FILE *ft, io_t io);
//~ extern int MBarmsol2(double *x,  mvbarms Prec);
//~ extern int condestMBArms(mvbarms armspre, double *y, FILE *fp );
//~ extern int cleanMBARMS(mvbarms ArmsPre);
//~ extern int setupilutpvector( iilutptr vbbmat, int nset, int *nBB);
//~ extern int computLUFp(vbsptr F, ilutptr lu, vbsptr FF);
//~ extern int computLUFinpartitionp( vbbsptr Bb, vbbsptr Ff, vbbsptr FF, iilutptr luu,io_t io, FILE *flog);
//~ extern int computschurpartitionp(vbp4ptr vbmat, vbsptr B, vbsptr C, double *droptol, int *lfil, vbsptr schur, int *nBB, int nset);
//~ extern int ilutpB( vbsptr schur, double *droptol, int *lfil, vbilutptr vbmat);
//~ extern void iilutptrsolver(iilutptr ilusch, double *y);
//~ extern int nnz_iilutptr(iilutptr lu );
//~ extern int cleanIILUT( iilutptr vbbmat );
//~ extern void qqsortm(int *ja,int *jab, BData *ma, int left, int right);
//~ extern int vb2csr2d(vbsptr vbmat, cs2dptr mat);
//~ extern int outputcs2dptr( cs2dptr mat, char *filename, int onebase);
//~ extern int csr2d2vb(cs2dptr mat, vbsptr vbmat, int n, int nc, int *bsz, int *bszc);
//~ extern int cleanCS2D(cs2dptr amat);
//~ extern int computLUFsparse(vbsptr F, vbiluptr lu, vbsptr FF);
int KeyComp( const void *vfst, const void *vsnd );
#endif 
