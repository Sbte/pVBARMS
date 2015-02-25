/**
 * @file   parms_vec.h
 * @author Zhongze Li
 * @date   Tue Oct 17 12:01:25 2006
 * 
 * @brief  Functions related to the parms_Vec object. 
 * 
 * 
 */

#ifndef _PARMS_VECTOR_H_
#define _PARMS_VECTOR_H_

#include "parms_sys.h"
#include "parms_map.h"
#include "parms_viewer.h"
#include "parms_comm.h"

PARMS_CXX_BEGIN

/** 
 * Return the 2-norm of the vector.
 * 
 * @param self  A vector object.    
 * @param value The 2-norm returned.
 * 
 * @return 0 on success.
 */
extern int parms_VecGetNorm2(FLOAT *self, REAL *value, parms_Map map);
extern int parms_VecGetNorm2_b(FLOAT *self, REAL *value, parms_Map map);

/** 
 * Scale a vector.
 *
 * All components of vector self on the local processor times scalar. 
 * \f$self = scalar \times self\f$.
 *  
 * @param self   A vector object.
 * @param scalar A scalar.      
 * 
 * @return 0 on success.
 */



extern int parms_VecScale(FLOAT *self, FLOAT scalar, parms_Map map);
extern int parms_VecScale_b(FLOAT *self, FLOAT scalar, parms_Map map);

/** 
 * Perform \f$self := scalar \times x + self\f$.
 * 
 * @param self   A vector object.      
 * @param x      Another vector object.
 * @param scalar A scalar.
 * 
 * @return 0 on success.
 */
extern int parms_VecAXPY(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map);
extern int parms_VecAXPY_b(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map);

/** 
 * Perform \f$self = scalar \times self + x\f$.
 * 
 * @param self    A vector object.      
 * @param x       Another vector object.
 * @param scalar  A scalar.
 * 
 * @return 0 on success.
 */
extern int parms_VecAYPX(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map);
extern int parms_VecAYPX_b(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map);

/** 
 * Perform the (global) inner product of two vectors.
 * 
 *  value = self x^{T}. If self
 * 
 *
 * @param self   A vector object.                  
 * @param x 	 Another vector object.            
 * @param value  The inner product returned.       
 * 
 * @return 0 on success.
 */
extern int parms_VecDOT(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map);
extern int parms_VecDOT_b(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map);
/** 
 * Perform the (global) inner product of two vectors.
 * 
 * If self and x are real vectors, value = self x^{T}. If self
 * and x are complex vectors, value = self \overline{x}^{T}.
 *
 * @param self   A vector object.                  
 * @param x 	 Another vector object.            
 * @param value  The inner product returned.       
 * 
 * @return 0 on success.
 */
int parms_VecDOTC(FLOAT *self, FLOAT *x, REAL *value, parms_Map is);
int parms_VecDOTC_b(FLOAT *self, FLOAT *x, REAL *value, parms_Map is);

/** 
 * Permute the vector object self.
 *
 * If the parms_Vec object and the parms_Mat object are created based
 * on  the same parms_Map object. Once matrix object is setup,
 * variables on each processor are divided into two cateries: internal
 * unknowns and interface unknowns. The parms_Vec should be permuted
 * accordingly. The user needn't call self function directly.
 * 
 * @param self A vector object.
 * 
 * @return 0 on success.
 */
extern int parms_VecPerm(FLOAT *self, parms_Map map);
extern int parms_VecPerm_b(FLOAT *self, parms_Map map);
/** 
 * Inverse permutation of the vector object self.
 * The user needn't call this function directly.
 * 
 * @param self A vector object.
 * 
 * @return 0 on success.
 */


extern int parms_VecInvPerm(FLOAT *self, parms_Map map);
extern int parms_VecInvPerm_b(FLOAT *self, parms_Map map);


extern int parms_VecSetValues_b(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map);
			       


PARMS_CXX_END

#endif 
