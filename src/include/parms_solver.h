/**
 * @file   parms_solver.h
 * @author Zhongze Li
 * @date   Tue Oct 17 11:58:22 2006
 * 
 * @brief  Functions related to the Krylov subspace methods. Only
 *         FGMRES is supported. 
 * 
 */

#ifndef _PARMS_SOLVER_H_
#define _PARMS_SOLVER_H_

#include "parms_vec.h"
#include "parms_mat.h"
#include "parms_viewer.h"
#include "parms_pc.h"
#include "parms_operator.h"

PARMS_CXX_BEGIN

typedef struct parms_Solver_ *parms_Solver;


extern int parms_SolverApply_b(parms_Solver self, FLOAT *x, FLOAT *y); 


/** 
 * Set the type of the solver.
 *
 * Only FGMRES solver is available in the package.
 * 
 * @param self   A parms_Solver object.       
 * @param stype  The type of Krylov subspace. 
 *               - SOLFGMRES
 *               - SOLDGMRES
 *               
 * @return 0 on success.
 */
extern int parms_SolverSetType(parms_Solver self, SOLVERTYPE stype);
extern int parms_SolverSetType_b(parms_Solver self, SOLVERTYPE stype);

/** 
 * Create a parms_Solver object.
 * 
 * @param self A pointer to the parms_Solver object created.
 * @param A    The matrix of the linear system.   
 * @param pc   The preconditioner.  
 * 
 * @return 0 on success.
 */
extern int parms_SolverCreate(parms_Solver *self, parms_Mat A, parms_PC pc);

extern int parms_SolverFree_b(parms_Solver *self);

/** 
 * Set parameter for the solver.
 * 
 * Set the maximum iteration counts, the restart size of GMRES, and
 * the convergence tolerance.
 * 
 * @param self   A parms_Solver object.
 * @param ptype  The type of parameter.
 *               - MAXITS maximum iteration counts.
 *               - KSIZE  restart size of GMRES.
 *               - DTOL   converence tolerance.
 *               - NEIG   number of eigenvectors.
 * @param param  Parameters for the solver.
 */
extern void parms_SolverSetParam(parms_Solver self, PARAMTYPE ptype,
				 char *param);
   
/** 
 * Get the iteration counts.
 * 
 * @param self A parms_Solver object.
 * 
 * @return The iteration counts.
 */
extern int parms_SolverGetIts(parms_Solver self);




PARMS_CXX_END

#endif 
