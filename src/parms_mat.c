/*--------------------------------------------------------------------
  parms_MatVec          : perform the multiplication of
                            matrix-vector product. 
  parms_MatCreate         : create a matrix object.
  parms_MatExtend         : extend submatrix by including equations
                            correspond to the immediate neighbouring
                            variables. 
  parms_MatFree           : free the memory for a matrix object.
  parms_MatFreeSubMat     : free the memory for the submatrix object.
  parms_MatGetCommHandler : get the communication handler.
  parms_MatGetDiag        : get the diagonal part of the local matrix. 
  parms_MatGetSubMat      : get the local matrix.
  parms_MatMVPY           : perform \f$z = \alpha \times self \times x
                                           + beta \times y\f$. 
  parms_MatSetCommType    : set the communication type.
  parms_MatSetValues      : insert/add values to the matrix object.
  parms_MatSetup          : set up a matrix object.
  parms_MatVecOffDiag     : perform the multiplication of the
                            off-diagonal matrix and the external vars. 
  parms_MatView           : dump a matrix object.

  A code fragment for using matrix functions:

  parms_Mat mat;
  parms_Map map;
  parms_PC  pc;

  ...

  // create a matrix object.
  parms_MatCreate(&mat, map);

  // insert/add values into the matrix
  parms_MatSetValues(mat, ...);

  // setup the matrix. 
  // this function divides local variables into two categories:
  // interior variables and interface variables. interior variables
  // listed first. when the matrix-vector product is performed, the
  // vector is permuted automatically if the vector is created
  parms_MatSetup(mat);

  // matrix-vector product
  parms_MatVec(mat, vec, y);

  // create a preconditioner object.
  parms_PCCreate(&pc, mat);

  $Id: parms_mat.c,v 1.4 2006-12-18 22:23:40 zzli Exp $
  ------------------------------------------------------------------*/
#include "include/parms_mat_impl.h"
#include "DDPQ/protos.h"//new

/*
int parms_MatCreate_vcsr(parms_Mat self);
int parms_MatCreate_dvcsr(parms_Mat self);
int parms_MatFree_dvcsr(parms_Mat *self);
int parms_MatView_vcsr(parms_Mat self, parms_Viewer v);
int parms_MatView_dvcsr(parms_Mat self, parms_Viewer v);
*/



/** 
 * Set the communication type.
 *
 * Set the communication style across processors.
 * communication style:
 *  -P2P       point-to-point (data copied to/from auxilliary buffers).
 *  -DERIVED   derived datatype.
 *
 * @param self  A matrix object.
 * @param ctype Communication style:
 *              - P2P     point-to-point (data copied to/from
 *                        auxilliary buffers).
 *              - DERIVED derived datatype.
 *              
 * @return 0 on success.
 */
int parms_MatSetCommType(parms_Mat self, COMMTYPE ctype)
{
  return self->ops->setcommtype(self, ctype);
}

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
int parms_MatCreate(parms_Mat *self, parms_Map map)
{
  parms_Mat  new_mat;

  PARMS_NEW0((new_mat));
  new_mat->ref = 1;
  PARMS_NEW0((new_mat)->ops);
  new_mat->is               = map;
  map->ref++;
  new_mat->type             = MAT_NULL;
  new_mat->ilutype          =PCILU0;     /*-- default local precon is ilu0. See parms_pc.c ---*/
  new_mat->issetup          = false;
  new_mat->isperm           = map->isperm;
  new_mat->isalloc          = false;
  new_mat->isreset          = false; 
  new_mat->isassembled        = false;   

  new_mat->m = parms_MapGetLocalSize(map);
  new_mat->n = new_mat->m;
  new_mat->M = parms_MapGetGlobalSize(map);
  new_mat->N = new_mat->M;
  new_mat->isserial = map->isserial;
  *self = new_mat;
  //--------------------------------------------------------------
  new_mat->isin_b_schur = false;/*turn it on when it is in block version schur routines*/


  return 0;
}

/** 
 * Perform \f$y = self \times x\f$.
 * 
 * @param self A parms_Mat object.      
 * @param x    A vector.      
 * @param y    Another vector object.
 * 
 * @return 0 on success.
 */
int parms_MatVec(parms_Mat self, FLOAT *x, FLOAT *y)
{
  return self->ops->apply(self, x, y);
}


/** 
 * Get the diagonal part of the local matrix. 
 * 
 * @param self A parms_Mat object.
 * @param mat  The diagonal part of the local matrix.    
 * 
 * @return 0 on success.
 */
int parms_MatGetDiag(parms_Mat self, void **mat)
{
  return self->ops->getdiag(self, mat);
}

/** 
 * Perform \f$z = alpha*self*x + beta*y\f$.
 * 
 * @param self   A matrix object.           
 * @param alpha  A scalar.                  
 * @param x      A vector object.           
 * @param beta   A scalar.                  
 * @param y      A vector object.           
 * @param z      A vector stores the result.
 * 
 * @return 0 on success.
 */
int parms_MatMVPY(parms_Mat self, FLOAT alpha, FLOAT *x, FLOAT beta,
		  FLOAT *y, FLOAT *z)
{
  return self->ops->mvpy(self, alpha, x, beta, y, z);
}

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
 *  function performs 
 *  \f[
 *    y[pos..n] = M_{ext} \times x_{ext}
 *  \f]
 *
 * @param self A matrix object.                                      
 * @param x    A vector object.                                      
 * @param y    A vector object.                                      
 * @param pos  The offset of x from the beginning of the local
 *             vector. 
 * 
 * @return 0 on success.
 */
int parms_MatVecOffDiag(parms_Mat self, FLOAT *x, FLOAT *y, int pos)
{
  return self->ops->mvoffd(self, x, y, pos);
}

/** 
 * Free the memory for the submatrix.
 * 
 * @param self A parms_Mat object.             
 * @param mat  The submatrix to be freed.      
 * 
 * @return 0 on success.
 */
int parms_MatFreeSubMat(parms_Mat self, void *mat)
{
  return self->ops->matfree(self, mat);
}

/** 
 * Extend submatrix by including equations correspond to the
 * immediate neighbouring variables.
 * 
 * @param self     A matrix object.                                         
 * @param handler  A communication handler.                                 
 * @param start    The beginning location of mat in the local matrix.       
 * @param mat      The submatrix to be extended.        
 * @param n 	 The size of extended matrix returned.
 * @param ext_mat  The extended matrix created.            
 * 
 * @return 0 on success.
 */
int parms_MatExtend(parms_Mat self, parms_Comm handler, int start,
		    void *mat, int *n, void **ext_mat)
{

  return self->ops->extend(self, handler, start, mat, n, ext_mat);
}

/** 
 * Get the local matrix. 
 * 
 * @param self A matrix object.                               
 * @param mat  The submatrix returned in a specific format.   
 * 
 * @return 0 on success.
 */
int parms_MatGetSubMat(parms_Mat self, void **mat)
{
  return self->ops->getlmat(self, mat);
}

/** 
 * Get the communication handler.
 * 
 * @param self    A matrix object.                     
 * @param handler The communication handler returned.  
 * 
 * @return 0 on success.
 */
int parms_MatGetCommHandler(parms_Mat self, parms_Comm *handler)
{
  return self->ops->gethandler(self, handler);
}

//----------------------------------------------------------------------------------------dividing line ----------------------------------------------------------
/** 
 * Insert/add values to the parms_Mat object self.
 * 
 * @param self    A parms_Mat object.                                
 * @param m 	The number of rows inserted.      
 * @param im 	An array of row global indices.                                                 
 * @param ia 	An array of pointer to the beginning of each row in
 *                array ja.    
 * @param ja      An array of column global indices.     
 * @param values  An array of values.                    
 * @param mode 	Insert value mode:                     
 * 		      - INSERT  insert values to parms_Mat self.
 * 		      - ADD     add values to parms_Mat self.
 *
 * NOTE: New entries are always inserted first, so mode does not
 *       matter if this is a new entry. Subsequent calls will either 
 *       replace (mode = INSERT) or add (mode = ADD) to the existing 
 *       entries in a particular position 
 * 		  
 * @return 0 on success.
 */
int parms_MatSetValues_b(parms_Mat self, int m, int *im, int *ia,
               int *ja, FLOAT **values, INSERTMODE mode)  //add mode is still missing, each entry of matrix point to the long array BData *values
{
  parms_Map  is;
  parms_bvcsr aux_data;
  BOOL isserial, isalloc;
  int i, j, k, offset, start, end, rindex, cindex;//row and column index
  int space, numinf, index, lrindex, lcindex, *rja;
  int *perm, *rp, *cp, *odp,inc, lsize, size, *nB, *lnB = NULL;//nB is global dimension of each block, lnB is local
  BData *ra;
  BOOL found;

  isalloc  = self->isalloc;
  is       = self->is;
  offset    = is->start;
  lsize    = parms_MapGetLocalSize(is);
  isserial = is->isserial;
  nB = is->nB;//for diagonal block size of local matrix

  if (isalloc) {
    aux_data = self->b_aux_data;     
  }
  else {
    self->isalloc = true;
    PARMS_NEW(aux_data);
    aux_data->n = lsize;
    PARMS_NEWARRAY(aux_data->space,   aux_data->n);
    PARMS_NEWARRAY0(aux_data->nzcount, aux_data->n);
    PARMS_NEWARRAY(aux_data->ja,      aux_data->n);
    PARMS_NEWARRAY(aux_data->ba,      aux_data->n);
    PARMS_NEWARRAY(aux_data->bsz,      aux_data->n+1);//vbsptr only
    //bsz[0] = 0;
    aux_data->bsz[0] = 0;
    aux_data->D = NULL;
    aux_data->bszc = NULL;
    aux_data->lsize = 0;



    for (i = 0; i < aux_data->n; i++) {
      aux_data->space[i] = 30;//why 30??
      PARMS_NEWARRAY(aux_data->ja[i], aux_data->space[i]); 
      PARMS_NEWARRAY(aux_data->ba[i], aux_data->space[i]);      
    }
    self->b_aux_data = aux_data;
    /* create table for tracking off-diagonal column count */
    parms_TableCreate(&self->odtable, NULL, lsize);
  }

  if (isserial) {
    for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      ra     = aux_data->ba[rindex];
      rja    = aux_data->ja[rindex];
      start  = ia[i] - offset;
      end    = ia[i+1] - offset;
      for (j = start; j < end; j++) {
	cindex = ja[j] - offset;
	space = aux_data->space[rindex];
	found = false;
	for (k = 0; k < aux_data->nzcount[rindex]; k++) {
	  if (rja[k] == cindex) {
	    found = true;
	    index = k;
	    break;
	  }
	}
	if (found) {
	  if (mode == INSERT) {
	    ra[index] = values[j];
	  }
	  else if (mode == ADD) {
	    //ra[index] += values[j];//modify it later,missing 
	  }
	}
	else { /* insert a new entry */
	  if (space == aux_data->nzcount[rindex]) {
	    /* reallocate memory for holding new entry */
	    space += 10;
	    PARMS_RESIZE(aux_data->ba[rindex], space);
	    PARMS_RESIZE(aux_data->ja[rindex], space);
	    aux_data->space[rindex] = space;	    
	    ra  = aux_data->ba[rindex];
	    rja = aux_data->ja[rindex];
	  }
	  rja[aux_data->nzcount[rindex]]  = cindex;
	  ra[aux_data->nzcount[rindex]++] = values[j];
	}
      }
    //output_intvector("nBinparms.coo",nB,0, m);getchar();
    aux_data->bsz[i+1] = aux_data->bsz[i] + nB[i];
    }
  }
  else {
    PARMS_NEWARRAY0(lnB,lsize);//allocate memory to array of local dimension of each block, only for parallel case, vbsptr only
    perm    = is->perm;
    numinf  = is->ninf;
/* Check if matrix is to be re-used with same non-zero pattern */
    if(self->isreset && (self->resetpattern == SAME_NONZERO_STRUCTURE))
    {
      for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      rp = parms_MapGlobalToLocal(is, rindex);
      /* this row resides the local processor */
      if (rp != NULL && *rp < lsize) {

	/* calculate the local row index */
	lrindex = perm[*rp];
	ra      = aux_data->ba[lrindex];
	rja     = aux_data->ja[lrindex];
        lnB[lrindex] = nB[rindex];//block and parallel case only, find the block size of the local row from the global row, vbsptr only
	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	for (j = start; j < end; j++) {
	  cindex = ja[j] - offset;
	  /* the local processor contains this column index */   
	  cp = parms_MapGlobalToLocal(is, cindex);
	  if (cp != NULL && *cp < lsize) {
	    /* local column index */
	    lcindex = perm[*cp];

	    found = false;
	    for (k = 0; k < aux_data->nzcount[lrindex]; k++) {
	      if (rja[k] == lcindex) {
	        found = true;
	        index = k;
	        break;
	      }
	     }
	     if (found) {
	       if (mode == INSERT) {
	         ra[index] = values[j];
	       }
	       else if (mode == ADD) {
	         //ra[index] += values[j];
	       }
	     }
	  }
	  else {
	    if(cp != NULL){
	      lcindex = *cp;
	    
	      found = false;
	      for (k = 0; k < aux_data->nzcount[lrindex]; k++) {
	        if (rja[k] == lcindex) {
	          found = true;
	          index = k;
	          break;
	        }
	       }
	       if (found) {
	         if (mode == INSERT) {
	           ra[index] = values[j];
	         }
	         else if (mode == ADD) {
	           //ra[index] += values[j];
	         }
	       }
	     }
	  }
	}
     }//if (rp != NULL && *rp < lsize)
     }//for (i = 0; i < m; i++)
     for (i = 0; i < lsize; i++) {//vbsptr
	 aux_data->bsz[i+1] = aux_data->bsz[i] + lnB[i];//vbsptr
     }//vbsptr
   }
   else{        
    for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
      rp = parms_MapGlobalToLocal(is, rindex);
      /* this row resides the local processor */
      if (rp != NULL && *rp < lsize) {
	/* calculate the local row index */
	lrindex = *rp;
	ra      = aux_data->ba[lrindex];
	rja     = aux_data->ja[lrindex];
        lnB[lrindex] = nB[rindex];//block and parallel case only, find the block size of the local row from the global row
	start = ia[i]-offset;
	end   = ia[i+1]-offset;
	for (j = start; j < end; j++) {
	  cindex = ja[j] - offset;
	  /* the local processor contains this column index */   
	  cp = parms_MapGlobalToLocal(is, cindex);
	  if (cp != NULL && *cp < lsize) {
	    /* local column index */
	    lcindex = *cp;
	  }
	  else {
	    if (perm[lrindex] == -1) { /* not marked yet */
	      perm[lrindex] = lsize-1-numinf;
	      numinf++;
	    }
	    lcindex = -cindex-1;
	  }
	  space   = aux_data->space[lrindex];
	  found = false;
	  for (k = 0; k < aux_data->nzcount[lrindex]; k++) {
	    if (rja[k] == lcindex) {
	      found = true;
	      index = k;
	      break;
	    }
	  }
	  if (found) {
	    if (mode == INSERT) {
	      ra[index] = values[j];
	    }
	    else if (mode == ADD) {
	      //ra[index] += values[j];
	    }
	  }
	  else { /* insert the new entry */
	    if (space == aux_data->nzcount[lrindex]) {
	      /* reallocate memory for holding new entry */
	      space += 10;
	      PARMS_RESIZE(aux_data->ba[lrindex], space);
	      PARMS_RESIZE(aux_data->ja[lrindex], space);
	      aux_data->space[lrindex] = space;
	      ra  = aux_data->ba[lrindex];
	      rja = aux_data->ja[lrindex];
	    }
	    if (cp == NULL) {
	      size = parms_TableGetSize(is->table);
	      parms_TablePut(is->table, cindex, size);
           /* Initialize column count */
	      inc = 0;
	      parms_TablePut(self->odtable,cindex,inc);
	    }
	    else if((cp !=NULL) && (lcindex < 0)){
	      /* update off-diagonal column count */	    
	       odp = parms_TableGet(self->odtable, cindex);
             inc = *odp + 1;
	       parms_TablePut(self->odtable,cindex,inc);
	    }
	    
	    rja[aux_data->nzcount[lrindex]]  = lcindex;
	    ra[aux_data->nzcount[lrindex]++] = values[j];
	   }
	  }
       }
      }
     is->ninf = numinf;
    }
    for (i = 0; i < lsize; i++) {//vbsptr
	aux_data->bsz[i+1] = aux_data->bsz[i] + lnB[i];//vbsptr
    }//vbsptr
  }
  is->llsize = aux_data->bsz[lsize];
  //is->lbsz_before = aux_data->bsz;
  PARMS_NEWARRAY(is->lbsz_before,lsize+1);//allocate memory to array of local dimension of each block, only for parallel case, vbsptr only
  PARMS_MEMCPY(is->lbsz_before, aux_data->bsz, lsize+1);
  if (lnB != NULL)
     { 
	free(lnB); 
	lnB = NULL;
     }
  //free(nB); nB = NULL;
  return 0;
}



/** 
 * Set up parms_Mat object self.
 *
 * This is the most important function for the parms_Mat object. This
 * function combines the function bdry and setup in the old version 
 * of pARMS. The function sets up the data structure needed by the
 * distributed matrix-vector multiplication, divides the variables on
 * the local processors into two categories: interior and interface
 * variables.
 * 
 * @param self A parms_Mat object. 
 * 
 * @return 0 on success.
 */
int parms_Mat_B_Setup(parms_Mat self)//block version
{      
  /* free and create new table for offdiagonal column count */
  if(self->odtable)
    parms_TableFree(&self->odtable);

  if ((!self->issetup) && (self->type == MAT_NULL)) {
    self->type = MAT_VCSR;
    if(self->isserial){
	parms_MatCreate_b_vcsr(self);
    }
    else{
      if(!self->isreset) 
	  parms_MatCreate_b_dvcsr(self);
    }
  }
  //printf("before b_dvcsr_setup\n");getchar();
  return self->ops->setup(self);
}




/** 
 * Free the parms_Mat object pointed to by self.
 * 
 * @param self A pointer to a parms_Mat object.
 * 
 * @return 0 on success.
 */
int parms_MatFree_b(parms_Mat *self) 
{
  int         m, i, nnz;

  (*self)->ref--;
  if ((*self)->ref == 0 ) {
    if(!(*self)->isserial)
	parms_MatFree_b_dvcsr(self);
    parms_MapFree(&(*self)->is);
    if ((*self)->isalloc) {
      m = (*self)->m;
      for (i = 0; i < m; i++) {
	  nnz = (*self)->b_aux_data->nzcount[i];
	  if (nnz) {
	    PARMS_FREE((*self)->b_aux_data->ja[i]);
	    PARMS_FREE((*self)->b_aux_data->ba[i]);
	  }
      }
//bsz and D
      PARMS_FREE((*self)->b_aux_data->ba);
      PARMS_FREE((*self)->b_aux_data->ja);
      PARMS_FREE((*self)->b_aux_data->nzcount);
      PARMS_FREE((*self)->b_aux_data->bsz);//new
      if ((*self)->b_aux_data->bszc) 
	PARMS_FREE((*self)->b_aux_data->bszc);//new
      if ((*self)->b_aux_data->D) {
        PARMS_FREE((*self)->b_aux_data->D);
      }
      PARMS_FREE((*self)->b_aux_data);
      PARMS_FREE((*self)->ops);
      PARMS_FREE(*self);
    }
  }
  return 0;
}
