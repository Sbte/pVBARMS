#ifndef __VBLOCK_HEADER_H__
#define __VBLOCK_HEADER_H__

#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#include "parms_mat_impl.h"
#include "parms_opt_impl.h"

#define MAX_BLOCK_SIZE   1000

/* FORTRAN style vblock format, compatible for many FORTRAN routines */
#define DATA(a,row,i,j)  (a[(j)*(row)+(i)])

/* the dimension of ith Block */
#define B_DIM(bs,i)      (bs[i+1]-bs[i])

#if 0
typedef struct SparRow {
    /*---------------------------------------------
| C-style CSR format - used internally
| for all matrices in CSR format 
|---------------------------------------------*/
    int n;
    int *nnzrow;  /* length of each row */
    int **ja;     /* pointer-to-pointer to store column indices  */
    double **ma;  /* pointer-to-pointer to store nonzero entries */
} SparMat, *csptr;
#endif 




typedef struct parms_vcsr SparMat;
typedef parms_vcsr csptr;

//typedef struct parms_bvcsr VBSparMat;
typedef struct Spa2Dmt {
    /*---------------------------------------------
| sparse 2D pointer to matrix,
| only store the nonzero columns. 
| specially for F in Barms
|---------------------------------------------*/
    int n;	/* row dimension of matrix */
    int nc;	/* column dimension of matrix */
    int nzcount;  /* length of each row */
    int *ja;     /* pointer to store column indices, sorted*/
    int *jab;	/* pointer to store block column indices, sorted*/
    FLOAT **ma;  /* pointer to store full dense column pointer */
} Spar2DMat, *cs2dptr;

typedef struct ILUfac {
    int n;
    csptr L;      /* L part elements                            */
    FLOAT *D;    /* diagonal elements                          */
    csptr U;      /* U part elements                            */
    int *work;    /* working buffer */
} ILUSpar, LDUmat, *iluptr;

typedef struct VBBSpaFmt {
    int n;        /* the block dimension of the independent sets*/
    int maxsize;
    int *bsz;     /* the row/col of the first element of each coarse block*/
    vbsptr *vba;   /* pointer-to-pointer to store coarse blocks */
} VBBSparMat, *vbbsptr;

typedef struct VBILUfac {
    int n;
    int *bsz;     /* the row/col of the first element of each   */
    /* diagonal block                             */
    BData *D;     /* diagonal blocks                            */
    vbsptr L;     /* L part blocks                              */
    vbsptr U;     /* U part blocks                              */
    int *work;    /* working buffer                             */
    BData bf;     /* buffer of a temp block                     */
    int DiagOpt;  /* Option for diagonal inversion/solutiob     */
    /* opt =  1 -->> call luinv                   */
    /* opt == 2 -->> block inverted call dgemv    */
} VBILUSpar, *vbiluptr; 

typedef struct VBBILUfac {
    int n;        /* the block dimension of the independent sets*/
    int *bsz;     /* the row/col of the first element of each coarse block*/
    vbiluptr *vlu;   /* pointer-to-pointer to store coarse blocks */
} VBBILUSpar, *vbbiluptr;

typedef struct PerMat4 *p4ptr;
typedef struct PerMat4 {
    /*------------------------------------------------------------
| struct for storing the block LU factorization 
| contains all the block factors except the 
| data related to the last block. 
| n       = size of current block
| symperm = whether or not permutations are symmetric.
|           used only in cleanP4..
| nB      = size of B-block
| L, U    = ILU factors of B-block
| F, E    = sparse matrices in (1,2) and (2,1) 
|           parts of matrix. 
| perm    = (symmetric) permutation used in factorization
|           comes from the independent set ordering
| rperm   = unsymmetric permutation (rows) not used in this
|           version -- but left here for compatibility..
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|----------------------------------------------------------*/ 
    int n;
    int nB;
    int symperm;
    /*   LU factors  */
    csptr L;
    csptr U;
    /* E, F blocks   */
    csptr E;
    csptr F;
    int *rperm;       /* row permutation         */
    int *perm;        /* col. permutation        */
    double *D1 ;      /* diagonal scaling row    */
    double *D2 ;      /* diagonal scaling columns*/
    FLOAT *wk;       /* work array              */
    /* pointer to next and previous struct         */
    p4ptr prev;
    p4ptr next;
} Per4Mat;

/* -------------------------------------------------------------------*/
typedef struct ILUTfac *ilutptr;
typedef struct ILUTfac {
    /*------------------------------------------------------------
| struct for storing data related to the last schur complement 
| we need to store the C matrix associated with the last block
| and the ILUT factorization of the related Schur complement.
| 
| n       = size of C block = size of Schur complement
| C       = C block of last level matrix. 
| L, U    = ILU factors of last schur complement. 
|
| meth[4] = parameters for defining variants in factorization 
|           - see function readin for details
| rperm    = row permutation used for very nonsymmetric matrices 
|            [such as bottleneck transversal] -- NOT IN THIS VERSION
| perm2     = unsymmetric permutation (columns) - used primarily
|           for the ILUTP version of ILUT/.. 
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|-----------------------------------------------------------*/
    int n;
    /*-------------------- C matrix of last block */
    csptr C;
    /* LU factorization       */
    csptr L;
    csptr U;
    /*--------------------  related to variants and methods */
    /*    int meth[4];   */
    int *rperm;   /* row single-sinded permutation */
    int *perm;    /* column perm .                */
    int *perm2;   /* column permutation coming from pivoting in ILU */
    double *D1;
    double *D2;
    FLOAT *wk;
} IluSpar;

typedef struct IILUTfac {//for coarse point ilutp struct
    int n;        /* the block dimension of the independent sets*/
    int *bsz;     /* the row/col of the first element of each coarse block*/
    ilutptr *llu;   /* pointer-to-pointer to store coarse blocks */
} IIluSpar, *iilutptr;

typedef struct VBILUTfac *vbilutptr;
typedef struct VBILUTfac {
    /*------------------------------------------------------------
| struct for storing data related to the last schur complement 
| we need to store the C matrix associated with the last block
| and the ILUT factorization of the related Schur complement.
| 
| n       = block dimension of C block = block dimension of Schur complement
| C       = C block of last level matrix. 
| L, U    = ILU factors of last schur complement. 
|
| meth[4] = parameters for defining variants in factorization 
|           - see function readin for details
| rperm    = row permutation used for very nonsymmetric matrices 
|            [such as bottleneck transversal] -- NOT IN THIS VERSION
| perm2     = unsymmetric permutation (columns) - used primarily
|           for the ILUTP version of ILUT/.. 
| D1, D2  = diagonal matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|-----------------------------------------------------------*/
    int n;
    int *bsz;
    /*-------------------- C matrix of last block */
    struct VBSpaFmt *C;
    /* LU factorization       */
    vbiluptr lu;
    ilutptr plu;
    //struct VBSpaFmt *L;
    //struct VBSpaFmt *U;
    //BData *D;
    /*--------------------  related to variants and methods */
    /*    int meth[4];   */
    int *rperm;   /* row single-sinded permutation */
    int *perm;    /* column perm .                */
    //int *perm2;   /* column permutation coming from pivoting in ILU */
    double *D1;
    double *D2;
    FLOAT *wk;
    double *dd1;//specially for if ( nB == 0 || nC == 0 ) case
    double *dd2;//specially for if ( nB == 0 || nC == 0 ) case
    //BData bf;     /* buffer of a temp block                     */
    //int DiagOpt;  /* Option for diagonal inversion/solutiob     */
    /* opt =  1 -->> call luinv                   */
    /* opt == 2 -->> block inverted call dgemv    */
} VBIluSpar;

typedef struct VBPerMat4 *vbp4ptr;
typedef struct VBPerMat4 {
    /*------------------------------------------------------------
| struct for storing the block LU factorization 
| contains all the block factors except the 
| data related to the last block. 
| n       = block dimension of current block
| symperm = whether or not permutations are symmetric.
|           used only in cleanP4..
| nB      = block dimension of B-block
| L, U    = ILU factors of B-block
| F, E    = sparse matrices in (1,2) and (2,1) 
|           parts of matrix. 
| perm    = (symmetric) permutation used in factorization
|           comes from the independent set ordering
| rperm   = unsymmetric permutation (rows) not used in this
|           version -- but left here for compatibility..
| D1, D2  = diagonal block matrices (left, right) used for scaling
|           if scaling option is turned on. Note that the 
|           method works by scaling the whole matrix first
|           (at any level) before anything else is done. 
| wk     = a work vector of length n needed for various tasks
|            [reduces number of calls to malloc]           
|----------------------------------------------------------
|              permuted and sorted matrices for each level of the block 
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
|            plus a few other things. See LIB/heads.h for details.*/
    int n;
    int *bsz; /* the row/col of the first element of each   */
    /* diagonal block                             */
    int nB;
    int symperm;
    /*   LU factors and D, block matrix  */
    /*BData *D ;
  struct VBSpaFmt *L;
  struct VBSpaFmt *U;*/
    vbiluptr lu;
    iilutptr plu;
    /* E, F blocks, block matrix   */
    struct VBSpaFmt *E;
    struct VBSpaFmt *F;
    int *rperm;       /* row permutation         */
    int *perm;        /* col. permutation        */
    double *D1 ;      /* block diagonal scaling row    */
    double *D2 ;      /* block diagonal scaling columns*/
    FLOAT *wk;       /* work array              */
    //BData bf;     /* buffer of a temp block                     */
    //int DiagOpt;  /* Option for diagonal inversion/solutiob     */
    /* opt =  1 -->> call luinv                   */
    /* opt == 2 -->> block inverted call dgemv    */
    /* pointer to next and previous struct         */
    vbp4ptr prev;
    vbp4ptr next;
} VBPer4Mat;

typedef struct parms_arms_data {
    int n;			/* dimension of matrix */
    int nlev;			/* number of levels */
    ilutptr ilus;
    p4ptr levmat;
    int ipar[18];
    double pgfpar[2];
    int schur_start;
    int ind;
    int nnz_mat;
    int nnz_prec;
} *parms_arms_data;

typedef struct parms_barms_data {
    int n;			/* dimension of matrix */
    int nlev;			/* number of levels */
    vbilutptr ilus;
    vbp4ptr levmat;
    vbilutptr ilusp;	/* ILU for last level, point ilutp only */
    int ipar[18];
    double pgfpar[2];
    int schur_start;//to store the point index of schur_start
    //int schur_start_p;//to store the point index of schur_start//changed
    int ind;
    int nnz_mat;
    int nnz_prec;
} *parms_barms_data;

typedef struct parms_arms_data armsMat;
typedef parms_arms_data arms;
typedef parms_barms_data vbarms;

typedef struct __CompressType
{
    int grp;   /* -1: begin new group, >=0: belong to grp-th row */
    int count; /* block size, valid only if grp = -1 */
} CompressType;

typedef struct __KeyType
{
    int var;   /* row number */
    int key;   /* hash value */
} KeyType;

#endif  /* __VBLOCK_HEADER_H__ */
