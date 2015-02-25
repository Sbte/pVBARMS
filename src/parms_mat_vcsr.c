#include "include/parms_opt_impl.h"
#include "include/parms_mat_impl.h"
#include "DDPQ/protos.h"//new


static int MatVec_b_vcsr(parms_Mat self, FLOAT *x, FLOAT *y)
{
  int          lsize, i, j, *ja, length, dim, sz, nBs, nBsj, col, *bsz, inc = 1;
  parms_bvcsr  matrix;
  parms_Map    is;
  BData        *ba;
  FLOAT       one = 1.0;
  
  
  /* extract diagonal and off-diagonal matrix */
  matrix   = self->b_aux_data;
  is       = self->is;
  lsize    = parms_MapGetLocalSize(is);
  //outputvbmatpa(matrix,"vbmatinmatvec.coo",1);//int outputvbmat( vbsptr vbmat, char *filename, int onebase)

  // perform local matrix vector product 
  bsz = matrix->bsz;
  for (i = 0; i < lsize; i++) {
    nBs = bsz[i];
    dim = B_DIM(bsz,i);

    for( j = 0; j < dim; j++ ) 
      y[nBs+j] = 0.0;

    length = matrix->nzcount[i];
    ja = matrix->ja[i];
    ba = matrix->ba[i];

    for (j = 0; j < length; j++) {
      col = ja[j];
      nBsj = bsz[col];
      sz = B_DIM(bsz,col);
      GGEMV ("n", dim, sz,one, ba[j],dim,&x[nBsj],inc,one,&y[nBs],inc);
    }
  }

  return 0;
}

static int MatMVPY_b_vcsr(parms_Mat self, FLOAT alpha, FLOAT *x,
			FLOAT beta, FLOAT *y, FLOAT *z) //not done yet
{
  //~ int          lsize, i, j, index, *pj, nnz;
  //~ parms_vcsr   matrix;
  //~ parms_Map    is;
  //~ FLOAT        *pa, s;
  //~ 
  //~ /* extract diagonal and off-diagonal matrix */
  //~ matrix   = self->aux_data;
  //~ is       = self->is;
  //~ lsize    = is->lsize;
//~ 
  //~ for (i = 0; i < lsize; i++) {
    //~ s = beta * y[i];
    //~ nnz = matrix->nnzrow[i];
    //~ pj  = matrix->pj[i];
    //~ pa  = matrix->pa[i];
    //~ for (j = 0; j < nnz; j++) {
      //~ index = pj[j];
      //~ s    += alpha * pa[j] * x[index];
    //~ }
    //~ z[i] = s;
  //~ }
  //~ 
  return 0;
}

static int MatSetup_b_vcsr(parms_Mat self)
{
  /* free member space in aux_data */
  //printf("in MatSetup\n");getchar();
  if(self->b_aux_data->space)
    PARMS_FREE(self->b_aux_data->space);
  self->issetup = true;
  self->is->ninf_send = 0;
  self->is->schur_start =  self->is->lsize;
  return 0;
}

static int MatGetDiag_b_vcsr(parms_Mat self, void **mat)//why dont use *csptr intead of void
{
  if (self->ilutype == PCVBARMS || self->ilutype == PCVBARMSOLD) {
    parms_bvcsr diag, diag_mat;
    int i, j, nnz, dimR, dimC, col, bnz;

    diag_mat = self->b_aux_data;
    PARMS_NEW(diag);
    diag->n = diag_mat->n;
    diag->bszc = NULL;
    diag->lsize = 0;
    diag->D = NULL;
    PARMS_NEWARRAY(diag->bsz, diag->n+1);
    diag->bsz[0] = 0;

    PARMS_NEWARRAY(diag->nzcount, diag->n);
    PARMS_NEWARRAY(diag->ja,     diag->n);
    PARMS_NEWARRAY(diag->ba,     diag->n);
    for (i = 0; i < diag->n; i++) {
      nnz = diag->nzcount[i] = diag_mat->nzcount[i];
      diag->bsz[i+1] = diag_mat->bsz[i+1];

      if (nnz) {
        dimR = B_DIM(diag_mat->bsz,i);
        PARMS_NEWARRAY(diag->ja[i], nnz);
        PARMS_NEWARRAY(diag->ba[i], nnz);
        PARMS_MEMCPY(diag->ja[i], diag_mat->ja[i], nnz);
        for (j = 0; j < nnz; j++){ 
            col = diag->ja[i][j];
            dimC = B_DIM(diag_mat->bsz,col);
            bnz = dimC*dimR;
            if (bnz) {
                PARMS_NEWARRAY(diag->ba[i][j], bnz);
                PARMS_MEMCPY(diag->ba[i][j], diag_mat->ba[i][j], bnz);
            }
        }
      }
    }
    *mat = diag;
    return 0;
  }
  *mat = self->b_aux_data;
  return 0;
}

static struct parms_Mat_ops parms_matops_b_vcsr = {
  MatVec_b_vcsr,
  MatSetup_b_vcsr,
  0,
  MatMVPY_b_vcsr,
  MatGetDiag_b_vcsr,
  MatGetDiag_b_vcsr,
  0,
  0,
  0,
  0
};
//~ 
//~ int parms_MatView_vcsr(parms_Mat self, parms_Viewer v)
//~ {
  //~ int        i, j, lsize, pid, length, *rja;
  //~ FLOAT      *ra;
  //~ parms_vcsr matrix;
  //~ parms_Map  is;
  //~ FILE       *fp;
//~ 
  //~ is      = self->is;
  //~ lsize   = parms_MapGetLocalSize(is);
  //~ pid     = is->pid;
  //~ parms_ViewerGetFP(v, &fp);
//~ 
  //~ fprintf(fp, "There is one processor available\n");
  //~ fprintf(fp, "The format of output local equations is:\n");
  //~ fprintf(fp, "(pid,local_row_index)=>(pid,global_row_index)\n");
  //~ fprintf(fp, "(pid,local_row_index, local_column_index) = value\n"); 
  //~ fprintf(fp, "\n");
//~ 
  //~ for (i = 0; i < lsize; i++) {
    //~ fprintf(fp, "(%d,%d) => (%d, %d)\n", pid, i, pid, i);
  //~ }
  //~ 
  //~ matrix = self->aux_data;
  //~ fprintf(fp, "Local matrix on processor %d is\n", pid);
//~ #if defined(DBL_CMPLX)  
  //~ for (i = 0; i < lsize; i++) {
    //~ length = matrix->nnzrow[i];
    //~ rja    = matrix->pj[i];
    //~ ra     = matrix->pa[i];
    //~ for (j = 0; j < length; j++) {
      //~ fprintf(fp, "(%d,%d,%d) = (%f, %f)  ", pid, i, rja[j], creal(ra[j]), cimag(ra[j]));
    //~ }
    //~ fprintf(fp, "\n");
  //~ }
//~ #else
  //~ for (i = 0; i < lsize; i++) {
    //~ length = matrix->nnzrow[i];
    //~ rja    = matrix->pj[i];
    //~ ra     = matrix->pa[i];
    //~ for (j = 0; j < length; j++) {
      //~ fprintf(fp, "(%d,%d,%d) = %f  ", pid, i, rja[j], ra[j]);
    //~ }
    //~ fprintf(fp, "\n");
  //~ }
//~ #endif
  //~ parms_ViewerStoreFP(v, fp);
  //~ return 0;
//~ }

//~ int parms_MatViewCOO_vcsr(parms_Mat self, parms_Viewer v)
//~ {
  //~ int        i, j, lsize, pid, length, *rja;
  //~ FLOAT      *ra;
  //~ parms_vcsr matrix;
  //~ parms_Map  is;
  //~ FILE       *fp;
//~ 
  //~ is      = self->is;
  //~ lsize   = parms_MapGetLocalSize(is);
  //~ pid     = is->pid;
  //~ parms_ViewerGetFP(v, &fp);
//~ 
  //~ fprintf(fp, "There is one processor available\n");//
  //~ fprintf(fp, "The format of output local equations is:\n");//
  //~ fprintf(fp, "(local_row_index local_column_index value\n"); //
  //~ fprintf(fp, "\n");//
//~ /*
  //~ for (i = 0; i < lsize; i++) {
    //~ fprintf(fp, "(%d,%d) => (%d, %d)\n", pid, i, pid, i);
  //~ }
//~ */  
  //~ matrix = self->aux_data;
  //~ fprintf(fp, "Local diagonal matrix on processor %d is\n", pid);//
  //~ fprintf(fp, "%d %d %d\n", lsize, lsize, nnzCS(matrix) );getchar();
//~ #if defined(DBL_CMPLX)  
  //~ for (i = 0; i < lsize; i++) {
    //~ length = matrix->nnzrow[i];
    //~ rja    = matrix->pj[i];
    //~ ra     = matrix->pa[i];
    //~ for (j = 0; j < length; j++) {
      //~ fprintf(fp, "%d  %d  %f  %f  \n", i, rja[j], creal(ra[j]), cimag(ra[j]));
    //~ }
    //~ fprintf(fp, "\n");
  //~ }
//~ #else
  //~ for (i = 0; i < lsize; i++) {
    //~ length = matrix->nnzrow[i];
    //~ rja    = matrix->pj[i];
    //~ ra     = matrix->pa[i];
    //~ for (j = 0; j < length; j++) {
      //~ fprintf(fp, "%d  %d  %f  \n", i, rja[j], ra[j]);
    //~ }
    //~ fprintf(fp, "\n");
  //~ }
//~ #endif
  //~ parms_ViewerStoreFP(v, fp);
  //~ return 0;
//~ }


int parms_MatCreate_b_vcsr(parms_Mat self)
{
  //printf("self->ops=%d",self->ops);getchar();
  PARMS_MEMCPY(self->ops, &parms_matops_b_vcsr, 1);
  return 0;
}
