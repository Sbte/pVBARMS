#ifndef _AUX_HEADER_INCLUDED_H
#define _AUX_HEADER_INCLUDED_H
#define  MAX_LINE  100
#define MAX_MAT 100

#include "parms.h"
#include <byteswap.h>//biswap_32
#include "../../src/DDPQ/protos.h"//new


#define BUFFLEN 200
#define STR2U(p, p1)				\
  if (1) {					\
    p = buf;					\
    while (' ' == *p) {				\
      p++;					\
    }						\
    p1 = malloc((strlen(p)+1)*sizeof(*p1));	\
    p2 = p1;					\
    while (' ' != *p) {				\
      *p2 = toupper(*p);			\
      ++p;					\
      ++p2;					\
    }						\
    *p2 = '\0';					\
  } else




typedef struct fprm_ {
  PCTYPE     pctype;
  PCILUTYPE  pcilutype;
  double     pgfpar[2];
  double     tolind;
  double     droptol[7];
  int        ipar[18];//ipar[2] is the bsize
  int        lfil[8];
  double     eps;/* For block version local solver(vbarms and vbilut) only, angle tolerance of init-block subroutine */
  int        cosine;
} *fprm;


typedef enum SYMMTYPE {SYMM, SKEW_SYMM, HERMITION} SYMMTYPE;

#if defined(FORTRAN_CAPS)
#define fread_param_         FREAD_PARAM
#define fset_pc_params_      FSET_PC_PARAMS
#define fset_solver_params_  FSET_SOLVER_PARAMS
#define fprm_free_           FPRM_FREE
#define wreadmtc_            WREADMTC
#define readmtc_             READMTC
#define csrcsc_              CSRCSC
#define aplb_                APLB
#define dse_                 DSE
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define fread_param_         fread_param__
#define fset_pc_params_      fset_pc_params__
#define fset_solver_params_  fset_solver_params__
#define fprm_free_           fprm_free__
#define wreadmtc_            wreadmtc__
#define readmtc_             readmtc__
#define csrcsc_              csrcsc__
#define aplb_                aplb__
#define dse_                 dse__
#elif !defined(FORTRAN_UNDERSCORE)
#define fread_param_         fread_param
#define fset_pc_params_      fset_pc_params
#define fset_solver_params_  fset_solver_params
#define fprm_free_           fprm_free
#define wreadmtc_            wreadmtc
#define readmtc_             readmtc
#define csrcsc_              csrcsc
#define aplb_                aplb
#define dse_                 dse
#endif

extern int  read_param(char *fname, char mname[MAX_MAT][MAX_LINE], fprm prm);
extern void set_pc_params(parms_PC pc, fprm prm);
extern void set_solver_params(parms_Solver solver, fprm prm);
extern void fread_param_(char *fname, fprm *prm, char *matrix, int *matlen, int len);
extern void fset_pc_params_(parms_PC *pc, fprm *prm);
extern void fset_solver_params_(parms_Solver *solver, fprm *prm);
extern void fprm_free_(fprm *prm);
extern void wreadmtc_(int *nmax, int *nzmax, int *job, char *fname,int
		      *len, double *a, int *ja, int *ia, double *rhs,
		      int *nrhs, char *guesol, int *nrow, int *ncol,
		      int *nnz, char *title, char *key, char *type,
		      int *ierr); 
extern void readmtc_(int *nmax, int *nzmax, int *job, char *fname, double 
		     *a, int *ja, int *ia, double *rhs, int *nrhs, char 
		     *guesol, int *nrow, int *ncol,int *nnz, char *title, 
		     char *key, char *type, int *ierr);
extern void csrcsc_(int *n, int *job, int *ipos, double *a, int *ja,
		    int *ia, double *ao, int *jao, int *iao);  //      subroutine csrcsc (n,job,ipos,a,ja,ia,ao,jao,iao)
extern void aplb_(int *nrow, int *ncol, int *job, double *a, int *ja,
		  int *ia, double *b, int *jb, int *ib, double *c, int
		  *jc, int *ic, int *nnzmax, int *iw, int *ierr); 
extern void dse_(int *n, int *ja, int *ia, int *ndom, int *riord, int
		 *dom, int *idom, int *mask, int *jwk, int *link);

extern void dse_(int *n, int *ja, int *ia, int *ndom, int *riord, int
		 *dom, int *idom, int *mask, int *jwk, int *link);
extern void coocsr_(int *nrow, int *nnz, double *a, int *ir, int *jc, double *ao, int *jao, int *iao);  //      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao)
//-----------------------------------dividing line-------------------------------------------------------------------
extern void coo2csptr(int n, int nnz, FLOAT *a, int *ir, int *jc, csptr mat);
extern void set_pc_params_b(parms_PC pc, fprm prm);
extern void output_intvector(char *filename,int *v,int i0, int i1);
extern void output_csrmatrix(char *filename,int *ia, int *ja, double *a, int n);
extern void output_dblvector(char *filename, FLOAT *v, int i0, int i1);
extern int outputcsmat ( csptr mat, char *filename, int onebase);
extern int colunms2csptr(int n, int *ia, int *ja, FLOAT *a, csptr mat);
extern int csptr2colunms(csptr mat, int *n, int *ia, int *ja);
extern int outputvbmat1( vbsptr vbmat, char *filename, int onebase);
extern int outputvbmat( vbsptr vbmat, char *filename, int onebase);
extern int vbsptr2colunms(vbsptr mat, int *bnnz, int *ia, int *ja, BData *a);
extern int outputcsmatpa ( csptr mat, char *filename, int onebase);
extern int output_intvectorpa(char *filename,int *v,int i0, int i1);
extern int output_dblvectorpa(char *filename,FLOAT *v,int i0, int i1);
extern int outputvbmatpa ( vbsptr mat, char *filename, int onebase);
extern int bincols2csptr(int n, int *nnzptr, int *ja, FLOAT *a, csptr mat);
extern double bswap_double(double a);
extern double byteswap_double(double v);
extern int fillinuppertrg(int n, int nnz, FLOAT *a, int *ja, int *ia, FLOAT *b, int *jb, int *ib, SYMMTYPE sytype);
extern int fullmatize(int n, int nnz, FLOAT *a, int *ja, int *ia, FLOAT *b, int *jb, int*ib, char *type);
extern int local_read_bin_data_from_indices(FILE *binfile, long int M, long int N, int nz, int nnzptr[], int ja[], FLOAT val[], int *indices, int nindices, int *gia);
extern int local_read_bin_data_b(FILE *binfile, long int M, long int N, int nz, int nnzptr[], int ja[], FLOAT val[], int *idom, int *dom, int *perm, int *nB, int nBlock, int *gia);
extern int local_read_bin_data(FILE *binfile, long int M, long int N, int nz, int ia[], int ja[], FLOAT val[], int *idom, int *dom, int *gia);

//extern int parms_MatSetValues_b(parms_Mat self, int m, int *im, int *ia, int *ja, BData *values, INSERTMODE mode);  //add mode is still missing, each entry of matrix point to the long array BData *values
//-----------------------------------dividing line-------------------------------------------------------------------

/*---------- complex routines ------------*/
#if defined(DBL_CMPLX)
#if defined(FORTRAN_CAPS)
#define wreadmtc_            WREADMTC
#define readmtc_             READMTC
#define csrcsc_              CSRCSC
#define aplb_                APLB
#elif defined(FORTRAN_DOUBLE_UNDERSCORE)
#define wreadmtc_            wreadmtc__
#define readmtc_             readmtc__
#define csrcsc_              csrcsc__
#define aplb_                aplb__
#elif !defined(FORTRAN_UNDERSCORE)
#define wreadmtc_            wreadmtc
#define readmtc_             readmtc
#define csrcsc_              csrcsc
#define aplb_                aplb
#endif
extern void zwreadmtc_(int *nmax, int *nzmax, int *job, char *fname,int
		      *len, complex double *a, int *ja, int *ia, complex double *rhs,
		      int *nrhs, char *guesol, int *nrow, int *ncol,
		      int *nnz, char *title, char *key, char *type,
		      int *ierr); 
extern void zreadmtc_(int *nmax, int *nzmax, int *job, char *fname, complex double 
		     *a, int *ja, int *ia, complex double *rhs, int *nrhs, char 
		     *guesol, int *nrow, int *ncol,int *nnz, char *title, 
		     char *key, char *type, int *ierr);
extern void zcsrcsc_(int *n, int *job, int *ipos, complex double *a, int *ja,
		    int *ia, complex double *ao, int *jao, int *iao);  
extern void zaplb_(int *nrow, int *ncol, int *job, complex double *a, int *ja,
		  int *ia, complex double *b, int *jb, int *ib, complex double *c, int
		  *jc, int *ic, int *nnzmax, int *iw, int *ierr); 

extern void zcoocsr_(int *nrow, int *nnz, complex double *a, int *ir, int *jc,
         complex double *ao, int *jao, int *iao);
//      subroutine coocsr(nrow,nnz,a,ir,jc,ao,jao,iao);

#endif

#endif 
