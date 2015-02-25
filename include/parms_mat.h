/**
 * @file   parms_mat.h
 * @author Zhongze Li
 * @date   Tue Oct 17 10:05:30 2006
 * 
 * @brief  Functions related to the matrix computations.
 * 
 * 
 */

#ifndef _PARMS_MAT_H_
#define _PARMS_MAT_H_

#include "parms_sys.h"
#include "parms_vec.h"
#include "parms_viewer.h"
#include "parms_operator.h"
//#include "../src/DDPQ/protos.h"//new???

PARMS_CXX_BEGIN

typedef struct parms_Mat_ *parms_Mat;//only 

/** 
 * Create a parms_Mat object.
 *
 * Create a parms_Mat object based on data distribution layout map.
 * 
 * @param self  A pointer to the parms_Mat object created. 
 * @param map   A parms_Map object, which describes the data
 *              distribution among processors.  
 * 
 * @return 0 on success.
 */
extern int parms_MatCreate(parms_Mat *self, parms_Map map);



extern int parms_MatFree_b(parms_Mat *self);
/** 
 * Perform \f$y = self \times x\f$.
 * 
 * @param self A parms_Mat object.      
 * @param x    A Vector object.      
 * @param y    Another Vector object.
 * 
 * @return 0 on success.
 */
extern int parms_MatVec(parms_Mat self, FLOAT *x, FLOAT *y);

/** 
 * Set the communication type.
 *
 * Set the communication style across processors.
 * communication style:
 *  - P2P       point-to-point (data copied to/from auxilliary buffers).
 *  - DERIVED   derived datatype.
 *
 * @param self  A matrix object.
 * @param ctype Communication style:
 *              - P2P     point-to-point (data copied to/from
 *                        auxilliary buffers).
 *              - DERIVED derived datatype.
 *              
 * @return 0 on success.
 */
extern int parms_MatSetCommType(parms_Mat self, COMMTYPE ctype);

/** 
 * Get the diagonal part of the local matrix. 
 * 
 * @param self A parms_Mat object.
 * @param mat  The diagonal part of the local matrix.    
 * 
 * @return 0 on success.
 */
extern int parms_MatGetDiag(parms_Mat self, void **mat);


/** 
 * Perform \f$z = alpha*self*x + beta*y\f$.
 * 
 * @param self   A matrix object.           
 * @param alpha  A scalar.                  
 * @param x 	 A vector object.           
 * @param beta 	 A scalar.                  
 * @param y 	 A vector object.           
 * @param z 	 A vector stores the result.
 * 
 * @return 0 on success.
 */
extern int parms_MatMVPY(parms_Mat self, FLOAT alpha, FLOAT *x, FLOAT
			 beta, FLOAT *y, FLOAT *z); 

/** 
 * Perform the multiplication of the off-diagonal matrix and the
 * external vars. 
 *
 * The local matrix can be written as follows:
 * 
 *  \f[
 *  \left(
 *  \begin{array}{ccc}
 *    B   &   E & 0\\
 *    F   &   C & M_{ext}
 *  \end{array}
 *  \right)
 *  \f],
 *  where \f$\left(\begin{array}{cc}
 *    B  &   E \\
 *    F  &   C 
 *  \end{array}
 *  \right)\f$ corresponds to the variables on the local PE. This
 *  function 
 *  performs 
 *  \f[
 *    y[pos..n] = M_{ext} \times x_{ext}
 *  \f]
 *
 * @param self A matrix object.                                      
 * @param x    A vector object.                                      
 * @param y    A vector object.                                      
 * @param pos  The offset of x from the beginning of he local vector.
 * 
 * @return 0 on success.
 */
extern int parms_MatVecOffDiag(parms_Mat self, FLOAT *x, FLOAT *y, int
			       pos);

/** 
 * Get the communication handler.
 * 
 * @param self    A matrix object.                     
 * @param handler The communication handler returned.  
 * 
 * @return 0 on success.
 */
extern int parms_MatGetCommHandler(parms_Mat self, parms_Comm
				   *handler);

/** 
 * Free the memory for the submatrix.
 * 
 * @param self A parms_Mat object.             
 * @param mat  The submatrix to be freed.      
 * 
 * @return 0 on success.
 */
extern int parms_MatFreeSubMat(parms_Mat self, void *mat);

/** 
 * Get the local matrix. 
 * 
 * @param self A matrix object.                               
 * @param mat  The submatrix returned in a specific format.   
 * 
 * @return 0 on success.
 */
extern int parms_MatGetSubMat(parms_Mat self, void **mat);

/** 
 * Extend submatrix by including equations correspond to the
 * immediate neighbouring variables.
 * 
 * @param self     A matrix object.                                         
 * @param handler  A communication handler.                                 
 * @param start    The beginning location of mat in the local matrix.       
 * @param mat      The submatrix to be extended.        
 * @param n 	   The size of extended matrix returned.
 * @param ext_mat  The extended matrix created.            
 * 
 * @return 0 on success.
 */
extern int parms_MatExtend(parms_Mat self, parms_Comm handler, int
			   start, void *mat, int *n, void **ext_mat);

/*
 *
 * Fortran Wrapper Functions 
 *
*/

extern void parms_matvec_(parms_Mat *self, FLOAT *x, FLOAT *y, int *ierr);  

extern void parms_matcreate_(parms_Mat *self, parms_Map *map, int *ierr);

extern void parms_matfree_(parms_Mat *self, int *ierr); 

extern void parms_matmvpy_(parms_Mat *self, FLOAT *alpha, FLOAT *x, FLOAT
		    *beta, FLOAT *y, FLOAT *z, int *ierr);

extern void parms_matsetcommtype_(parms_Mat *self, COMMTYPE *ctype, int
			   *ierr);

extern void parms_matsetvalues_(parms_Mat *self, int *m, int *im, int *ia,
			 int *ja, FLOAT *values, INSERTMODE *mode, int
			 *ierr);
			 
extern void parms_matsetelementmatrix_(parms_Mat *self, int *m, int *im, int *ia,
			 int *ja, FLOAT *values, INSERTMODE *mode, int *ierr);

extern void parms_matassembleelementmatrix_(parms_Mat *self, int *ierr);	

extern void parms_matresetrowvalues_(parms_Mat *self, int *m, int *im, int *ia,
		       int *ja, FLOAT *values, int *ierr);	       

extern void parms_matreset_(parms_Mat *self, NNZSTRUCT *nonzerostructure, int *ierr);

extern void parms_matsetup_(parms_Mat *self, int *ierr);

extern void parms_matview_(parms_Mat *self, parms_Viewer *v, int *ierr);

extern void parms_matviewcoo_(parms_Mat *self, parms_Viewer *v, int *ierr);     

extern int parms_MatSetValues_b(parms_Mat self, int m, int *im, int *ia, int *ja, FLOAT **values, INSERTMODE mode);  //add mode is still missing, each entry of matrix point to the long array BData *values

//-----------------------------------------------------dividing line-----------------------------------------------    

extern int parms_Mat_B_Setup(parms_Mat self);
/*
 *
 * end Fortran Wrapper Functions 
 *
*/

PARMS_CXX_END

#endif 
