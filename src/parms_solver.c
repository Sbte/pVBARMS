/*--------------------------------------------------------------------
  parms_SolverApply    : solve the linear system of equation
  parms_SolverCreate   : create a parms_Solver object.
  parms_SolverFree     : free the memory for the parms_Solver object.
  parms_SolverGetIts   : get the iteration counts.
  parms_SolverGetMat   : get the matrix for the linear system solved.
  parms_SolverGetPC    : get the preconditioning matrix.
  parms_SolverSetParam : set the parameters for the solver (maxits,
                         tolerance, etc.)
  parms_SolverSetType  : set the type of the solver. FGMRES and DGMRES
                         are supplied in the package.
                         dgmres removed YS
  parms_SolverView     : dump the maximum iteration count and the
                         iteration counts

  A code fragment for using solver functions:

  parms_Solver solver;
  parms_Mat A;
  parms_PC  pc;

  // create a solver
  parms_SolverCreate(&solver, A, pc);

  // set the maximum number of iterations
  parms_SolverSetParam(solver, MAXITS, "300");
  // set the restart size of FGMRES or DGMRES
  parms_SolverSetParam(solver, KSIZE, "60");
  // set the convergence tolerance
  parms_SolverSetParam(solver, DTOL, "1.0e-6");
  // set the number of eigenvectors
  parms_SolverSetParam(solver, NEIG, "8");

  // solver the linear system
  parms_SolverApply(solver, rhs, x);

  // free the memory for solver
  parms_SolverFree(&solver);

  $Id: parms_solver.c,v 1.5 2006-12-01 20:44:20 zzli Exp $
 -------------------------------------------------------------------*/
#include <stdlib.h>
#include "parms_vec.h"
#include "parms_viewer.h"
#include "parms_mat_impl.h"
#include "parms_pc_impl.h"
#include "parms_solver_impl.h"


/** 
 * Create a parms_Solver object.
 *
 * @param self A pointer to the parms_Solver object created.
 * @param A    The matrix of the linear system.
 * @param pc   The preconditioner.
 *
 * @return 0 on success.
 */

int parms_SolverCreate(parms_Solver *self, parms_Mat A, parms_PC pc)
{
    parms_Solver new_solver;

    PARMS_NEW0((new_solver));
    new_solver->ref = 1;
    PARMS_NEW0((new_solver)->ops);
    new_solver->ops->apply    = 0;
    new_solver->ops->setksize = 0;
    new_solver->istypeset = false;
    new_solver->A = A;
    A->ref++;
    new_solver->pc = pc;
    pc->ref++;
    new_solver->maxits = 100;
    new_solver->tol    = 1.0e-6;
    *self = new_solver;
    return 0;
}





/** 
 * Set the type of the solver.
 *
 * Only FGMRES solver is available in the package.
 *
 * @param self   A parms_Solver object.
 * @param stype  The type of Krylov subspace.
 *               -SOLFGMRES
 *               -SOLDGMRES
 *
 * @return 0 on success.
 */

int parms_SolverSetType(parms_Solver self, SOLVERTYPE stype)
{
    if (self->istypeset && self->stype == stype) {
        return 0;
    }
    if (self->istypeset) {
        self->ops->solver_free(&self);
    }
    if(stype == SOLFGMRES)
        fgmres_create(self);
    else if(stype == SOLGMRES)
        gmres_create(self);
    else{
        printf("ERROR: Invalid choice of solver - (Check SOLVERTYPE for parms_SolverSetType(...) \n");
        PARMS_ABORT(17);
    }
    self->stype      = stype;
    self->istypeset  = true;
    return 0;
}

/** 
 * Set parameter for the solver.
 *
 * Set the maximum iteration counts, the restart size of GMRES, and
 * the convergence tolerance.
 *
 * @param self   A parms_Solver object.
 * @param ptype  The type of parameter.
 *               -MAXITS maximum iteration counts.
 *               -KSIZE  restart size of GMRES.
 *               -DTOL   converence tolerance.
 *               -NEIG   number of eigenvectors.
 * @param param  Parameters for the solver.
 */
void parms_SolverSetParam(parms_Solver self, PARAMTYPE paramtype, char
                          *param)
{

    if (self->istypeset == false) {
        parms_SolverSetType(self, SOLFGMRES);
    }
    if (paramtype == MAXITS) {
        self->maxits = atoi(param);
    }
    else if (paramtype == KSIZE) {
        self->ops->setksize(self, atoi(param));
    }
    else if (paramtype == DTOL) {
        self->tol = strtod(param, (char **)NULL);
    }
    else if (paramtype == NEIG) {
        self->ops->setneig(self, atoi(param));
    }
}


/** 
 * Get the iteration counts.
 *
 * @param self A parms_Solver object.
 *
 * @return The iteration counts.
 */
int parms_SolverGetIts(parms_Solver self)
{
    return self->its;
}






//------------------------------------------------------------------------------dividing line------------------------------------------------------------------------
/** 
 * Set the type of the solver.
 *
 * Only FGMRES solver is available in the package.
 *
 * @param self   A parms_Solver object.
 * @param stype  The type of Krylov subspace.
 *               -SOLFGMRES
 *               -SOLDGMRES
 *
 * @return 0 on success.
 */

int parms_SolverSetType_b(parms_Solver self, SOLVERTYPE stype)
{
    if (self->istypeset && self->stype == stype) {
        return 0;
    }
    if (self->istypeset) {
        self->ops->solver_free(&self);
    }
    if(stype == SOLFGMRES)
        fgmres_create_b(self);
    //else if(stype == SOLGMRES)
    //gmres_create(self);
    else{
        printf("ERROR: Invalid choice of solver - (Check SOLVERTYPE for parms_SolverSetType(...) \n");
        PARMS_ABORT(17);
    }
    self->stype      = stype;
    self->istypeset  = true;
    return 0;
}



/** 
 * Solve the equation \f$Ax = y\f$.
 *
 * @param self A parms_Solver object.
 * @param x    The solution vector.
 * @param y    The right-hand-side vector.
 *
 * @return 0 on success.
 */
int parms_SolverApply_b(parms_Solver self, FLOAT *y, FLOAT *x)
{
    //  if(self->istypeset == true){
    //      self->ops->solver_free(&self);
    //      self->istypeset = false;
    //  }
    if (self->istypeset == false) {
        /* Default solver - fgmres */
        parms_SolverSetType_b(self, SOLFGMRES);
    }

    return self->ops->apply(self, y, x);
}

/** 
 * Free the memory for the parms_Solver object.
 *
 * @param self A pointer to the parms_Solver object to be freed.
 *
 * @return 0 on success.
 */
int parms_SolverFree_b(parms_Solver *self)
{
    (*self)->ref--;
    if ((*self)->ref == 0 ) {
        parms_MatFree_b(&(*self)->A);
        parms_PCFree_b(&(*self)->pc);
        (*self)->ops->solver_free(self);
        PARMS_FREE((*self)->ops);
        PARMS_FREE(*self);
    }
    return 0;
}
