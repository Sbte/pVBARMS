#ifndef _PARMS_OPERATOR_IMPL_H_
#define _PARMS_OPERATOR_IMPL_H_

#include "parms_vec.h"
#include "parms_operator.h"
#include "parms_mat_impl.h"

typedef struct parms_Operator_ops{
 int  (*apply)(parms_Operator self, FLOAT *y, FLOAT *x);
 int  (*lsol)(parms_Operator self, FLOAT *y, FLOAT *x);
 int  (*invs)(parms_Operator self, FLOAT *y, FLOAT *x);
 int  (*ascend)(parms_Operator self, FLOAT *y, FLOAT *x);
 int  (*getssize)(parms_Operator self);
 void (*getnnz)(parms_Operator self, int *nnz_mat, int
		       *nnz_pc);
 int  (*operator_free)(parms_Operator *self);
 int  (*operator_view)(parms_Operator self, parms_Viewer v); 
} *parms_Operator_ops; 

struct parms_Operator_ {
  int ref;
  parms_Operator_ops ops;
  void *data;
  //BOOL islocalgc;//new, for judging if we use local graph compression or the global one.
  int *blockperm;//permutation array returned from init_block routine
};

/*-----external function protos--------*/
extern int parms_ilu0_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
extern int parms_iluk_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
extern int parms_ilut_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
extern int parms_barms_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
extern int parms_barms_vcsr_old(parms_Mat self, parms_FactParam param, void *mat, parms_Operator *op);
extern int parms_barms_vcsr_l(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
extern int parms_barms_vcsr_l_old(parms_Mat self, parms_FactParam param, void *mat, parms_Operator *op);
extern int parms_vbilut_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 
extern int parms_vbiluk_vcsr(parms_Mat self, parms_FactParam param, void *data, parms_Operator *op); 

/*----------End Protos------------------*/

#endif 
