static int MatVec_b_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y) 
{
 int lsize, i, j, *ja, length, dim, sz, nBs, nBsj, col, inc = 1, row, llsize;
 int *bsz, *bszc;
  parms_dvcsr data;
  parms_Comm  handler;
  parms_bvcsr b_diag_mat, b_offd_mat;
  parms_Map is;
  FLOAT *offsetptr;//, s;//*ba, 
  BData *ba;
  double one=1.0;

//printf("one = %20.16e\n", one);//getchar(); 
   //output_dblvectorpa("xinmatvec",x,0, self->is->llsize);
  parms_VecPerm_b(x, self->is);
  parms_VecPerm_b(y, self->is);
//printf("x = %d\n", x);//getchar(); 
  // extract diagonal and off-diagonal matrix 
  data     = (parms_dvcsr)self->data;
  handler  = data->mvhandler;
  b_diag_mat = data->b_diag_mat;
  b_offd_mat = data->b_offd_mat;
  is       = self->is;
  lsize    = is->lsize;
  llsize   = is->llsize;

  parms_CommDataBegin_b(handler, x, 0);

  // perform local matrix vector product 
  bsz = b_diag_mat->bsz;
  for (i = 0; i < lsize; i++) {
    //s = 0.0;
    nBs = bsz[i];
    dim = B_DIM(bsz,i);
    for( j = 0; j < dim; j++ ) 
      y[nBs+j] = 0;
    //memset(&y[nBs], 0, dim*sizeof(FLOAT));
    length = b_diag_mat->nzcount[i];
    ja = b_diag_mat->ja[i];
    ba = b_diag_mat->ba[i];
    for (j = 0; j < length; j++) {
      col = ja[j];
      nBsj = bsz[col];
      //index = ja[j];
      sz = B_DIM(bsz,col);
      //-------------------- operation:  y = Block*x + y 
      DGEMV ("n", dim, sz,one, ba[j],dim,&x[nBsj],inc,one,&y[nBs],inc);//s    += ba[j] * x[index];
    }
    //y[i] = s;
  }
   //output_dblvectorpa("yafterlocalproduct",y,0,self->is->llsize );//handler->ptrvrecv_b[handler->nprecv]
  // off-diagonal part matrix-vector product 
  parms_CommDataEnd(handler);//block and point version are the same on this routine
  offsetptr = handler->buf_recv_b;// - llsize;// because in the b_offd_mat the stating column index is above lsize.

//printf("is->ninf = %d, b_offd_mat->n = %d, b_offd_mat->nc = %d, offsetptr = %d\n", is->ninf, b_offd_mat->n , b_offd_mat->nc, offsetptr);//getchar(); 
   //output_dblvectorpa("offsetptr",offsetptr,0, 3);//handler->ptrvrecv_b[handler->nprecv]

  //MPI_Barrier(MPI_COMM_WORLD);
  //bsz = b_offd_mat->bsz;
  bszc = b_offd_mat->bszc;

  for (i = 0; i < is->ninf; i++) {

    row = i+is->nint;
    nBs = bsz[row];
    dim = B_DIM(bsz,row);


    length = b_offd_mat->nzcount[i];
    ja     = b_offd_mat->ja[i];
    ba     = b_offd_mat->ba[i];
    //s      = y[i+is->nint];// the stating row index of offd_mat is 0, so we need to plus the number of internal nodes

    for (j = 0; j < length; j++){ 
      col = ja[j]-lsize;
      nBsj = bszc[col];
      //index = ja[j];
      sz = B_DIM(bszc,col);
      //-------------------- operation:  y = Block*x + y 
      DGEMV ("n", dim, sz,one, ba[j],dim,&offsetptr[nBsj],inc,one,&y[nBs],inc);//s    += ba[j] * x[index];

      //s   += ba[j] * offsetptr[ja[j]];
//printf("i = %d, j = %d, ja[j] = %d, length = %d \n", i, j , ja[j], length);getchar(); 
    }
    //y[i+is->nint] = s;
  }
   //output_dblvectorpa("yinmatvecbeforepermutation",y,0,self->is->llsize );//handler->ptrvrecv_b[handler->nprecv]
  parms_VecInvPerm_b(x, self->is);
  parms_VecInvPerm_b(y, self->is);
   //output_dblvectorpa("yinmatvec",y,0,self->is->llsize );//handler->ptrvrecv_b[handler->nprecv]
  return 0;
}









int lsize, i, j, index, *ja, length, llsize, dim, sz, nBs, nBsj, col, inc = 1;
  parms_dvcsr data;
  parms_Comm handler;
  parms_bvcsr b_diag_mat, b_offd_mat;
  parms_Map is;
  FLOAT *offsetptr, s;// *pa,
  BData *ba;

/* no need to permute z here, since it only stores the result */  
  parms_VecPerm_b(x, self->is);
  parms_VecPerm_b(y, self->is);

  /* extract diagonal and off-diagonal matrix */
  data     = (parms_dvcsr)self->data;
  handler  = data->mvhandler;
  b_diag_mat = data->b_diag_mat;
  b_offd_mat = data->b_offd_mat;
  is       = self->is;
  lsize    = is->lsize;
  llsize    = is->llsize;

  /* post receive and send actions */
  parms_CommDataBegin_b(handler, x, 0);  

  bsz = b_diag_mat->bsz;
  /* perform local matrix vector product */
  for (i = 0; i < lsize; i++) {
    s = beta * y[i];//GSCAL(lsize, scalar, self, i);
    length = b_diag_mat->nzcount[i];
    ja = b_diag_mat->ja[i];
    ba = b_diag_mat->ba[i];
    for (j = 0; j < length; j++) {
      index = ja[j];
      s    += alpha * ba[j] * x[index];
    }
    z[i] = s;
  }

  parms_CommDataEnd(handler);

  offsetptr = handler->buf_recv - lsize;
  for (i = 0; i < is->ninf; i++) {
    length = b_offd_mat->nzcount[i];
    ja     = b_offd_mat->ja[i];
    ba     = b_offd_mat->ba[i];
    s      = z[i+is->nint];
    for (j = 0; j < length; j++) {
      s   += alpha * ba[j] * offsetptr[ja[j]];
    }
    z[i+is->nint] = s;
  }

/* inverse permutation. z now has permutation of y, so has to be inverted as well */
  parms_VecInvPerm_b(x, self->is);
  parms_VecInvPerm_b(y, self->is);
  parms_VecInvPerm_b(z, self->is);



















/*--------------------------------------------------------------------
  parms_VecAXPY          : y = y + alpha x
  parms_VecAYPX          : y = scalar y + x
  parms_VecDOT           : inner product of two vectors.
  parms_VecDotArray      : inner products for the vector to an array
                           of vector object.

  parms_VecGetNorm2      : get the 2-norm of a vector.
  parms_VecPerm          : permute a vector object.

  parms_VecScale         : scale the components of a vector.
  A code fragment for using vector functions:

  parms_Vec vec;
  parms_Map map;

  // assume PE0 contains global variables (0, 2, 4, 6, 8), PE1
  // contains global variables (1, 3, 5, 7, 9). Both PE0 and PE1 call
  // the following parms_VecSetValues function. Each PE will pick and
  // insert variables which belong to the local PE according to the
  // partionting information stored in map.

  // all vector computation functions can be called
  parms_VecGetNorm2(vec, &value, map);
  parms_VecDot(vec, vec, &value, map);
  ...

  $Id: parms_vec.c,v 1.3 2007-05-10 14:11:21 zzli Exp $
  ------------------------------------------------------------------*/
#include "parms_sys.h"
#include "parms_viewer.h"
#include "parms_vec.h"
#include "parms_map_impl.h"
#if defined(__ICC)
#include <mathimf.h>
#else
#if defined(C99)
#include <tgmath.h>
#else
#include <math.h>
#endif 
#endif  

typedef struct ext_vec_data {
  /*! \var im The external variables in global indexing
   */
  int     *im;
 /*!  \var values The external variable contribution.
  */
  FLOAT    *values;
 /*!  \var n Number if external variable contribution.
  */
  int    n; 
 /*!  \var space Size of work space.
  */
  int    space;    
} *ext_vec_data;


/** 
 * Scale a vector.
 *
 * All components of parms_Vec object self times scalar. 
 * \f$self = scalar \times self\f$.
 *  
 * @param self   A vector object.
 * @param scalar A scalar.      
 * 
 * @return 0 on success.
 */
int parms_VecScale(FLOAT *self, FLOAT scalar, parms_Map map) 
{
  int lsize, i;

  lsize = parms_MapGetLocalSize(map);
#ifdef HAS_BLAS
  i = 1;
  GSCAL(lsize, scalar, self, i);
#else
  for (i = 0; i < lsize; i++) 
  {
    self[i] *= scalar;
  }
#endif 
  return 0;
}


/** 
 * Perform \f$self := scalar \times x + self\f$.
 * 
 * @param self   A vector object.      
 * @param x      Another vector object.
 * @param scalar A scalar.
 * 
 * @return 0 on success.
 */
int parms_VecAXPY(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map)   
{
  int i, lsize;

  lsize = parms_MapGetLocalSize(map);

#ifdef HAS_BLAS
  i = 1;
  GAXPY(lsize, scalar, x, i, self, i);
#else
  for (i = 0; i < lsize; i++) 
 {
    self[i] += scalar * x[i];
  }
#endif 

  return 0;
}

/** 
 * Perform \f$self = scalar \times self + x\f$.
 * 
 * @param self    A vector object.      
 * @param x 	  Another vector object.
 * @param scalar  A scalar.
 * 
 * @return 0 on success.
 */
int parms_VecAYPX(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map) 
{
  int lsize, i;
  lsize = parms_MapGetLocalSize(map);

#ifdef HAS_BLAS
  FLOAT  t;
  i = 1;
  t = 1.0;
  GSCAL(lsize, scalar, self, i);
  GAXPY(lsize, t, x, i, self, i);
#else

  for (i = 0; i < lsize; i++) {
    self[i] = scalar * self[i] + x[i];
  }
#endif 
  return 0;
}


/* Dot product on local PE 
*/
static int vec_LDot(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map)    
{
  int lsize, i;

  lsize = parms_MapGetLocalSize(map);
#ifdef HAS_BLAS

#if defined(DBL_CMPLX)
  i = 1;
  *value = GDOTU(lsize, self, i, x, i);
#else
  i = 1;
  *value = GDOT(lsize, self, i, x, i);
#endif

#else 

#if defined(DBL_CMPLX)
  FLOAT dot = 0.0 + 0.0*I;
  for (i = 0; i < lsize; i++)
  	dot  +=  self[i] * x[i];  
#else
  FLOAT dot = 0.0;  
  for (i = 0; i < lsize; i++) {
      dot  +=  self[i] * x[i];
#endif

  *value = dot;
#endif 

  return 0;
}

/* Dot product on local PE - uses complex conjugate
*/
static int vec_LDotc(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map)    
{
  int lsize, i;

  lsize = parms_MapGetLocalSize(map);

#ifdef HAS_BLAS

#if defined(DBL_CMPLX)
  i = 1;
  *value = GDOTC(lsize, self, i, x, i);
#else
  i = 1;
  *value = GDOT(lsize, self, i, x, i);
#endif

#else 

#if defined(DBL_CMPLX)
  FLOAT dot = 0.0 + 0.0*I;
  for (i = 0; i < lsize; i++)
  	dot  +=  self[i] * conj(x[i]);
#else
  FLOAT dot = 0.0;  
  for (i = 0; i < lsize; i++) {
      dot  +=  self[i] * x[i];
#endif

  *value = dot;
#endif 

  return 0;
}

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
int parms_VecDOT(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map is)
{
  MPI_Comm comm;

  comm = is->comm;

  if (is->isserial) {
    vec_LDot(self, x, value, is);
  }
  else {
    FLOAT ldot;
    vec_LDot(self, x, &ldot, is);
    
#if defined(DBL_CMPLX)
    MPI_Allreduce(&ldot, value, 1, MPI_CMPLX, MPI_CMPLX_SUM, comm);
#else
    MPI_Allreduce(&ldot, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif    
  }

  return 0;
}

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
int parms_VecDOTC(FLOAT *self, FLOAT *x, REAL *value, parms_Map is)
{
  MPI_Comm comm;

  comm = is->comm;
  
  if (is->isserial) {
#if defined(DBL_CMPLX)  
    FLOAT val;    
    vec_LDotc(self, x, &val, is);
    *value = creal(val);    
#else
    vec_LDot(self, x, value, is);
#endif   
  }
  else {
#if defined(DBL_CMPLX)  
    FLOAT ldot;
    REAL lval;
    vec_LDotc(self, x, &ldot, is);
    lval = creal(ldot);

    MPI_Allreduce(&lval, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#else
    FLOAT ldot;
    vec_LDot(self, x, &ldot, is);
    MPI_Allreduce(&ldot, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif    
  }

  return 0;
}

/** 
 * Return the 2-norm of the vector.
 * 
 * @param self  A vector object.    
 * @param value The 2-norm returned.
 * 
 * @return 0 on success.
 */
int parms_VecGetNorm2(FLOAT *self, REAL *value, parms_Map is)
{ 

  if (is->isserial) {
#ifdef HAS_BLAS
    int lsize = parms_MapGetLocalSize(is);
    int incr = 1;
    *value = GNRM2(lsize, self, incr);
#else
    FLOAT dot;
    vec_LDotc(self, self, &dot, is);
    *value = ABS_VALUE(dot);
    *value = sqrt(*value);
#endif 
  }
  else {
    REAL dot;
#if defined(DBL_CMPLX)    
    parms_VecDOTC(self, self, &dot, is);
#else
    parms_VecDOT(self, self, &dot, is);
#endif 
    *value = fabs(dot);
    *value = sqrt(*value);
  }
  return 0;
}

/** 
 * Perform the inner product between self and an array of vectors.  
 *
 * The pseudo code:
 *
 *  \f{verbatim}
 *  for (i = 0; i < n; i++) {
 *    result[i] = self * vecarray[i];
 *  }
 *  \f}
 *
 * @param self       A vector object.                              
 * @param n 	     The size of vecarray.                         
 * @param vecarray   An array of vector objects.                   
 * @param aux 	     An auxiliary array.                           
 * @param result     An array of size n to store inner products.   
 * 
 * @return 0 on success.
 */
int parms_VecDotArray(FLOAT *self, int n, FLOAT
			**vecarray, FLOAT *result, parms_Map is)
{
  int i;
  FLOAT *aux;
  MPI_Comm comm;

  comm = is->comm;
  if (is->isserial) {
    for (i = 0; i < n; i++)
      vec_LDotc(self, vecarray[i], &result[i], is);
  }
  else {
    PARMS_NEWARRAY(aux, n);
    for (i = 0; i < n; i++) 
      vec_LDotc(self, vecarray[i], &aux[i], is);
#if defined(DBL_CMPLX)
      MPI_Allreduce(aux, result, n, MPI_CMPLX, MPI_CMPLX_SUM, comm);
#else
      MPI_Allreduce(aux, result, n, MPI_DOUBLE, MPI_SUM, comm);
#endif    
    PARMS_FREE(aux); 
  }
  return 0;
}

int parms_VecPerm(FLOAT *self, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  FLOAT        *newvec;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm && (map->isvecperm == false)) { 
    perm = map->perm;
    PARMS_NEWARRAY(newvec, size);
    for (i = 0; i < size; i++) {
      k = perm[i];
      newvec[k] = self[i];
    }
    memcpy(self, newvec, size*sizeof(FLOAT));

    PARMS_FREE(newvec);
  }
  return 0;
}

int parms_VecPermAux(FLOAT *self, FLOAT *aux, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm) { 
    perm = map->perm;
    for (i = 0; i < size; i++) {
      k = perm[i];
      aux[k] = self[i];
    }
  }
  else
  {
    for (i = 0; i < size; i++) {
      aux[i] = self[i];
    }
   }     
  return 0;
}

int parms_VecInvPerm(FLOAT *self, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  FLOAT        *newvec;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm && (map->isvecperm == false)) { 
    perm = map->iperm;
    PARMS_NEWARRAY(newvec, size);
    for (i = 0; i < size; i++) {
      k = perm[i];
      newvec[k] = self[i];
    }
    memcpy(self, newvec, size*sizeof(FLOAT));

    PARMS_FREE(newvec);
  }
  return 0;
}

int parms_VecInvPermAux(FLOAT *self, FLOAT *aux, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  size = parms_MapGetLocalSize(map);
  if (map->isperm) { 
    perm = map->iperm;
    for (i = 0; i < size; i++) {
      k = perm[i];
      aux[k] = self[i];
    }
  }
  else
  {
    for (i = 0; i < size; i++) {
      aux[i] = self[i];
    }
   }  
  return 0;
}

/** 
 * Insert values to parms_Vec object self.
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables to be inserted.     
 * @param im     An array of global variable indices.     
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
int parms_VecSetValues(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map) 
{
  int lsize, rindex, index, lrindex, *rp, i;
  int offset = map->start;

  lsize = parms_MapGetLocalSize(map);

  if (!map->isserial) {
    if (map->isvecperm) {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((index = *rp) < lsize)) {
	  lrindex = map->perm[index];
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if (mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      }
    }
    else {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((lrindex = *rp) < lsize)) {
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if(mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      }
    }
  }
  else {
    for (i = 0; i < m; i++) {
      rindex = im[i] - offset;
//printf("map->isvecperm=%d\n",map->isvecperm);
      if (map->isvecperm) {
	rindex = map->perm[rindex];
      }
      if (mode == ADD) {
	self[rindex] += values[i];
      }
      else if(mode == INSERT) {
	self[rindex] = values[i];
      }
    }
  }
  return 0;
}

/** 
 * Insert values to parms_Vec object self. This assumes the vector 
 * values are being set element-by-element. A call to parms_VecAssembleElementVector
 * is required to complete the vector once all entries have been added.
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables to be inserted.     
 * @param im     An array of global variable indices.     
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
int parms_VecSetElementVector(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map) 
{
  int lsize, rindex, index, lrindex, *rp, i,pos,k;
  int offset = map->start, space;
  ext_vec_data vdata;

  lsize = parms_MapGetLocalSize(map);
  
  if(map->isserial) {
/* call vecsetvalues for serial version */
    return parms_VecSetValues(self, m, im, values, mode, map);  
  }  

/* Parallel call */
  if(map->isdatalloc)
  {
    vdata = (ext_vec_data)map->data;
  }
  else
  {
  /* Allocate some memory for vdata */
    map->isdatalloc = true;
    PARMS_NEW(vdata);
    vdata->space = lsize;
    PARMS_NEWARRAY(vdata->im, vdata->space);
    PARMS_NEWARRAY0(vdata->values, vdata->space);
    vdata->n = 0;
    map->data = vdata;
  }

  if (map->isvecperm) {
    for (i = 0; i < m; i++) {
     rindex = im[i] - offset;
     rp = parms_MapGlobalToLocal(map, rindex);
     if ((rp != NULL) && ((index = *rp) < lsize)) {
	  lrindex = map->perm[index];
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if (mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      else{
	/* external contribution - first check if there is enough memory*/
	  space = vdata->space;
        if(vdata->n == space){
	  /* reallocate memory for holding new entry */
	  space += 10;
        PARMS_RESIZE(vdata->im,   space);
        PARMS_RESIZE(vdata->values, space);
        vdata->space = space;
        }	  

        pos = 0;
        for(k = 0; k<vdata->n; k++) /* check if row has been added already */
        {
          if(vdata->im[k] == im[i])
              break;
        }
        if(k < vdata->n) // already added row
        {
          pos = k;
          vdata->values[pos] += values[i];
        }
        else          // new row contribution
        {
          pos = vdata->n;
          vdata->im[vdata->n] = im[i];
          vdata->values[vdata->n++] = values[i];
        }
	}
     }
    }
    else {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((lrindex = *rp) < lsize)) {
	  if (mode == ADD) {
	    self[lrindex] += values[i];
	  }
	  else if(mode == INSERT) {
	    self[lrindex] = values[i];
	  }
	}
      else{
	/* external contribution - first check if there is enough memory*/
	  space = vdata->space;
        if(vdata->n == space){
	  /* reallocate memory for holding new entry */
	  space += 10;
        PARMS_RESIZE(vdata->im,   space);
        PARMS_RESIZE(vdata->values, space);
        vdata->space = space;
        }	  

        pos = 0;
        for(k = 0; k<vdata->n; k++) /* check if row has been added already */
        {
          if(vdata->im[k] == im[i])
              break;
        }
        if(k < vdata->n) // already added row
        {
          pos = k;
          vdata->values[pos] += values[i];
        }
        else          // new row contribution
        {
          pos = vdata->n;
          vdata->im[pos] = im[i];
          vdata->values[pos] = values[i];
          vdata->n++;
        }
	}
    }
   }
 
  return 0;
}

/** 
 * Completes setting up values for the distributed vector
 *
 * @param self   A vector object.                          
 * @param map    A pARMS map object
 *               
 * @return 0 on success.
 */
int parms_VecAssembleElementVector(FLOAT *self, parms_Map map) 
{
  FLOAT *gvec;
  int i, k, gnodv,lsize, mypid, npro, lrindex, index;
  int *info, *disp, *pos, *gextvars;
  MPI_Comm comm;
  ext_vec_data vdata;
  
  lsize = map->lsize;
  mypid = map->pid;
  npro = map->npro;
  comm = map->comm;

  if(map->isdatalloc){
    vdata = (ext_vec_data)map->data;
/* Allocate some memory */
    PARMS_NEWARRAY(info, npro);
    PARMS_NEWARRAY(disp, npro+1);

/* collect contributions from other processors */
    MPI_Allgather(&vdata->n, 1, MPI_INT, info, 1, MPI_INT, comm);
    disp[0] = 0;
    for (i = 0; i < npro; i++) {
      disp[i+1] = disp[i] + info[i];
    }
  
    gnodv = disp[npro];
    PARMS_NEWARRAY(gextvars, gnodv);
    PARMS_NEWARRAY0(gvec, gnodv);
    
    MPI_Allgatherv(vdata->im, vdata->n, MPI_INT, gextvars,
		   info, disp, MPI_INT, comm);
#if defined(DBL_CMPLX)
    MPI_Allgatherv(vdata->values, vdata->n, MPI_CMPLX, gvec,
		   info, disp, MPI_CMPLX, comm);   
#else
    MPI_Allgatherv(vdata->values, vdata->n, MPI_DOUBLE, gvec,
		   info, disp, MPI_DOUBLE, comm);  
#endif    
    for(i=0; i<npro; i++)
    {
      if(mypid != i)
      {
        for(k = disp[i]; k<disp[i+1]; k++)
        {
           index = gextvars[k] - map->start;
           pos = parms_MapGlobalToLocal(map,index);
           if ((pos != NULL) && ((lrindex = *pos) < lsize)) {
             self[lrindex] += gvec[k];
           }
         }
       }
     }    
/* Free memory */
    PARMS_FREE(info);
    PARMS_FREE(disp);
    PARMS_FREE(gvec);
    PARMS_FREE(gextvars);
    PARMS_FREE(vdata->im);
    PARMS_FREE(vdata->values);
    PARMS_FREE(vdata);
    map->isdatalloc = false;
  } 
   
  return 0;
}

/** 
 * Gather distributed vector to a global array.
 * 
 * @param self The distributed vector.
 * @param ga A global vector.
 * 
 * @return 0 on success.
 */
int parms_VecGather(FLOAT *self, FLOAT *ga, parms_Map map)
{
  BOOL isserial;
  int i, j, incx, lsize, index, npro, pid;
  int *num_recv, *ind_array, *ind_snd, maxnum;
  FLOAT *wk;

  isserial = map->isserial;
  pid      = map->pid;
  npro     = map->npro;

  lsize = parms_MapGetLocalSize(map);
  if (isserial) {
      incx = 1;
      GCOPY(lsize, self, incx, ga, incx);
  }
  else {
    PARMS_NEWARRAY0(num_recv, npro);
    PARMS_NEWARRAY(ind_snd,   lsize);
    MPI_Allgather(&lsize, 1, MPI_INT, num_recv, 1, MPI_INT, map->comm);
    maxnum = 0;
    for (i = 0; i < npro; i++) {
      if (maxnum < num_recv[i]) {
	maxnum = num_recv[i];
      }
    }
    PARMS_NEWARRAY(wk, maxnum);
    PARMS_NEWARRAY(ind_array, maxnum);

    for (i = 0; i < lsize; i++) {
	index = map->lvars[i];
	ga[index] = self[i];
	ind_snd[i] = map->lvars[i];
    }

#if defined(DBL_CMPLX)    
    for (i = 0; i < npro; i++) {
      if (pid == i) {
	MPI_Bcast(ind_snd, lsize, MPI_INT, i, map->comm);
	MPI_Bcast(self, lsize, MPI_CMPLX, i, map->comm);
      }
      else {
	MPI_Bcast(ind_array, num_recv[i], MPI_INT, i, map->comm);
	MPI_Bcast(wk, num_recv[i], MPI_CMPLX, i, map->comm);
	for (j = 0; j < num_recv[i]; j++) {
	  index = ind_array[j];
	  ga[index] = wk[j];
	}
      }
    }
#else
    for (i = 0; i < npro; i++) {
      if (pid == i) {
	MPI_Bcast(ind_snd, lsize, MPI_INT, i, map->comm);
	MPI_Bcast(self, lsize, MPI_DOUBLE, i, map->comm);
      }
      else {
	MPI_Bcast(ind_array, num_recv[i], MPI_INT, i, map->comm);
	MPI_Bcast(wk, num_recv[i], MPI_DOUBLE, i, map->comm);
	for (j = 0; j < num_recv[i]; j++) {
	  index = ind_array[j];
	  ga[index] = wk[j];
	}
      }
    }
#endif    
    PARMS_FREE(ind_snd);
    PARMS_FREE(ind_array);
    PARMS_FREE(wk);
  }

  return 0;
}


//----------------------------------------------------------------------------------------------dividing line-----------------------------------------------------------
int parms_VecPerm_b(FLOAT *self, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  FLOAT        *newvec;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //size = parms_MapGetLocalSize(map);
  size = parms_MapGetLocallSize(map);

  if (map->isperm && (map->isvecperm == false)) { 
    //PARMS_NEWARRAY(pperm, size);

    perm = map->pperm;
    PARMS_NEWARRAY(newvec, size);
    for (i = 0; i < size; i++) {
      k = perm[i];
      newvec[k] = self[i];
    }
    memcpy(self, newvec, size*sizeof(FLOAT));

    PARMS_FREE(newvec);
    //PARMS_FREE(pperm);
    //free(pperm);pperm = NULL;
  }
  return 0;
}

int parms_VecInvPerm_b(FLOAT *self, parms_Map map)
{
  int          *perm;
  int          size, i, k, rank;
  FLOAT        *newvec;
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //size = parms_MapGetLocalSize(map);
  size = parms_MapGetLocallSize(map);

  if (map->isperm && (map->isvecperm == false)) { 
    //PARMS_NEWARRAY(pperm, size);

    perm = map->pperm;
    PARMS_NEWARRAY(newvec, size);
    for (i = 0; i < size; i++) {
      k = perm[i];
      newvec[i] = self[k];
    }
    memcpy(self, newvec, size*sizeof(FLOAT));

    PARMS_FREE(newvec);
    //PARMS_FREE(pperm);
    //free(pperm);pperm = NULL;
  }
  return 0;
}


/** 
 * Insert values to parms_Vec object self. (block version)
 * 
 *  A pseudo code from the global point of view:
 *
 *  \f{verbatim}
 *  for (i = 0; i < m; i++) {
 *    self[im[i]] = values[i]; 
 *  }
 *  \f}
 *  
 * @param self   A vector object.                          
 * @param m      The number of variables(blocks) to be inserted.     
 * @param im     An array of global variable indices.     
 * @param value  An array of values to be inserted to self.s 
 * @param mode   The style of set values:
 *               -ADD    add values to parms_Vec object self. 
 *               -INSERT assign values to parms_Vec object self.
 *               
 * @return 0 on success.
 */
int parms_VecSetValues_b(FLOAT *self, int m, int *im, FLOAT
			       *values, INSERTMODE mode, parms_Map map) //not strictly tested yet
{
  int lsize, rindex, index, lrindex, *rp, i, j;
  int offset = map->start;
  int *nB, *lbsz_before, *lbsz_after, *gbsz;

  nB = map->nB;
  lbsz_before = map->lbsz_before;
  lbsz_after = map->lbsz_after;
  PARMS_NEWARRAY(gbsz, m+1);//vbsptr only
  gbsz[0] = 0;
  //PARMS_NEWARRAY(lnB, map->llsize);//vbsptr only
//printf("m=%d\n",m);
  for (i = 0; i < m; i++) {
	gbsz[i+1] = gbsz[i]+nB[i];
  }

  lsize = parms_MapGetLocalSize(map);

  if (!map->isserial) {
    if (map->isvecperm) {
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;
	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((index = *rp) < lsize)) {
	  lrindex = map->perm[index];
	  if (mode == ADD) {
		for (j = 0; j < nB[rindex]; j++) {
		   self[lbsz_after[lrindex]+j] += values[gbsz[rindex]+j];
		}	    //self[lrindex] += values[i];
	  }
	  else if (mode == INSERT) {
	    	for (j = 0; j< nB[rindex]; j++){
		   self[lbsz_after[lrindex]+j] = values[gbsz[rindex]+j];
		}//self[lrindex] = values[i];
	  }
	}
      }
    }
    else {
	printf("in insert and unpermuted mode \n");
      for (i = 0; i < m; i++) {
       rindex = im[i] - offset;

	rp = parms_MapGlobalToLocal(map, rindex);
	if ((rp != NULL) && ((lrindex = *rp) < lsize)) {

	  //lnB[lrindex] = nB[rindex];

	  if (mode == ADD) {
		for (j = 0; j< nB[rindex]; j++){
		   self[lbsz_before[lrindex]+j] += values[gbsz[rindex]+j];
		}
	    //self[lrindex] += values[i];
	  }
	  else if(mode == INSERT) {
                

		for (j = 0; j< nB[rindex]; j++){
		   self[lbsz_before[lrindex]+j] = values[gbsz[rindex]+j];
		}
	    //self[lrindex] = values[i];
	  }
	}
      }
    }
  }
  /*else {//serial case is not tested yet, I need to figure out what is isvecperm for in this case.
    //permuten(map->lvars, fulperm, A->b_aux_data->bsz, nloc);//int permuten(int *iwork, int *iworkn, int *bsz, int n)//here insert the expansion of permutation
    for (i = 0; i < m; i++) {
      //bim->pim(point version)
      rindex = im[i] - offset;
      if (map->isvecperm) {
	rindex = map->perm[rindex];
      }
      if (mode == ADD) {
		for (j = 0; j< nB[rindex]; j++){
		   self[lbsz_before[lrindex]+j] += values[i+j];
		}//self[rindex] += values[i];
      }
      else if(mode == INSERT) {
		for (j = 0; j< nB[rindex]; j++){
		   self[lbsz_before[lrindex]+j] += values[i+j];
		}//self[rindex] = values[i];
      }
    }
  }*/
  PARMS_FREE(gbsz);
  return 0;
}


/* Dot product on local PE - uses complex conjugate
*/
static int vec_LDotc_b(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map)    
{
  int llsize, i;

  llsize = parms_MapGetLocallSize(map);

#ifdef HAS_BLAS

#if defined(DBL_CMPLX)
  i = 1;
  *value = GDOTC(llsize, self, i, x, i);
#else
  i = 1;
  *value = GDOT(llsize, self, i, x, i);
#endif

#else 

#if defined(DBL_CMPLX)
  FLOAT dot = 0.0 + 0.0*I;
  for (i = 0; i < llsize; i++)
  	dot  +=  self[i] * conj(x[i]);
#else
  FLOAT dot = 0.0;  
  for (i = 0; i < llsize; i++) {
      dot  +=  self[i] * x[i];
#endif

  *value = dot;
#endif 

  return 0;
}

/* Dot product on local PE 
*/
static int vec_LDot_b(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map map)    
{
  int llsize, i;

  llsize = parms_MapGetLocalSize(map);
#ifdef HAS_BLAS

#if defined(DBL_CMPLX)
  i = 1;
  *value = GDOTU(llsize, self, i, x, i);
#else

  printf("in vec_LDot_b\n");
  i = 1;
  *value = GDOT(llsize, self, i, x, i);
#endif

#else 

#if defined(DBL_CMPLX)
  FLOAT dot = 0.0 + 0.0*I;
  for (i = 0; i < llsize; i++)
  	dot  +=  self[i] * x[i];  
#else
  FLOAT dot = 0.0;  
  for (i = 0; i < llsize; i++) {
      dot  +=  self[i] * x[i];
#endif

  *value = dot;
#endif 

  return 0;
}




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
int parms_VecDOTC_b(FLOAT *self, FLOAT *x, REAL *value, parms_Map is)
{
  MPI_Comm comm;

  comm = is->comm;
  
  if (is->isserial) {
#if defined(DBL_CMPLX)  
    FLOAT val;    
    vec_LDotc_b(self, x, &val, is);
    *value = creal(val);    
#else
    vec_LDot_b(self, x, value, is);
#endif   
  }
  else {
#if defined(DBL_CMPLX)  
    FLOAT ldot;
    REAL lval;
    vec_LDotc_b(self, x, &ldot, is);
    lval = creal(ldot);

    MPI_Allreduce(&lval, value, 1, MPI_DOUBLE, MPI_SUM, comm);//int MPI_Allreduce ( void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm )

#else
    FLOAT ldot;
    vec_LDot_b(self, x, &ldot, is);
    MPI_Allreduce(&ldot, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif    
  }

  return 0;
}


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
int parms_VecDOT_b(FLOAT *self, FLOAT *x, FLOAT *value, parms_Map is)
{
  MPI_Comm comm;

  comm = is->comm;

  if (is->isserial) {
    vec_LDot_b(self, x, value, is);
  }
  else {
    FLOAT ldot;
    vec_LDot_b(self, x, &ldot, is);
    
#if defined(DBL_CMPLX)
    MPI_Allreduce(&ldot, value, 1, MPI_CMPLX, MPI_CMPLX_SUM, comm);
#else
    MPI_Allreduce(&ldot, value, 1, MPI_DOUBLE, MPI_SUM, comm);
#endif    
  }

  return 0;
}

/** 
 * Return the 2-norm of the vector.(block version)
 * 
 * @param self  A vector object.    
 * @param value The 2-norm returned.
 * 
 * @return 0 on success.
 */
int parms_VecGetNorm2_b(FLOAT *self, REAL *value, parms_Map is)
{ 

  if (is->isserial) {
#ifdef HAS_BLAS
    int llsize = parms_MapGetLocallSize(is);
    int incr = 1;
    *value = GNRM2(llsize, self, incr);
#else
    FLOAT dot;
    vec_LDotc_b(self, self, &dot, is);
    *value = ABS_VALUE(dot);
    *value = sqrt(*value);
#endif 
  }
  else {
    REAL dot;
#if defined(DBL_CMPLX)    
    parms_VecDOTC_b(self, self, &dot, is);
#else
    parms_VecDOT_b(self, self, &dot, is);
#endif 
    *value = fabs(dot);
    *value = sqrt(*value);
  }
  return 0;
}


/** 
 * Perform \f$self := scalar \times x + self\f$.
 * 
 * @param self   A vector object.      
 * @param x      Another vector object.
 * @param scalar A scalar.
 * 
 * @return 0 on success.
 */
int parms_VecAXPY_b(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map)   
{
  int i, lsize;

  lsize = parms_MapGetLocallSize(map);//should be point lsize

#ifdef HAS_BLAS
  i = 1;
  GAXPY(lsize, scalar, x, i, self, i);
#else
  for (i = 0; i < lsize; i++) 
 {
    self[i] += scalar * x[i];
  }
#endif 

  return 0;
}


/** 
 * Perform \f$self = scalar \times self + x\f$.(block version)
 * 
 * @param self    A vector object.      
 * @param x 	  Another vector object.
 * @param scalar  A scalar.
 * 
 * @return 0 on success.
 */
int parms_VecAYPX_b(FLOAT *self, FLOAT *x, FLOAT scalar, parms_Map map) 
{
  int lsize, i;
  lsize = parms_MapGetLocallSize(map);//should be point lsize instead of block lsize

#ifdef HAS_BLAS
  FLOAT  t;
  i = 1;
  t = 1.0;
  GSCAL(lsize, scalar, self, i);
  GAXPY(lsize, t, x, i, self, i);
#else

  for (i = 0; i < lsize; i++) {
    self[i] = scalar * self[i] + x[i];
  }
#endif 
  return 0;
}

/** 
 * Scale a vector.
 *
 * All components of parms_Vec object self times scalar. 
 * \f$self = scalar \times self\f$.
 *  
 * @param self   A vector object.
 * @param scalar A scalar.      
 * 
 * @return 0 on success.
 */
int parms_VecScale_b(FLOAT *self, FLOAT scalar, parms_Map map) 
{
  int lsize, i;

  lsize = parms_MapGetLocallSize(map);//should be point lsize

//printf("in parms_VecScale_b\n");

#ifdef HAS_BLAS
  i = 1;
  GSCAL(lsize, scalar, self, i);
#else
  for (i = 0; i < lsize; i++) 
  {
    self[i] *= scalar;
  }
#endif 
  return 0;
}


/** 
 * Perform the inner product between self and an array of vectors. (block version) not tested yet 
 *
 * The pseudo code:
 *
 *  \f{verbatim}
 *  for (i = 0; i < n; i++) {
 *    result[i] = self * vecarray[i];
 *  }
 *  \f}
 *
 * @param self       A vector object.                              
 * @param n 	     The size of vecarray.                         
 * @param vecarray   An array of vector objects.                   
 * @param aux 	     An auxiliary array.                           
 * @param result     An array of size n to store inner products.   
 * 
 * @return 0 on success.
 */
int parms_VecDotArray_b(FLOAT *self, int n, FLOAT
			**vecarray, FLOAT *result, parms_Map is)
{
  int i;
  FLOAT *aux;
  MPI_Comm comm;

  comm = is->comm;
  if (is->isserial) {
    for (i = 0; i < n; i++)
      vec_LDotc_b(self, vecarray[i], &result[i], is);
  }
  else {
    PARMS_NEWARRAY(aux, n);
    for (i = 0; i < n; i++) 
      vec_LDotc_b(self, vecarray[i], &aux[i], is);
#if defined(DBL_CMPLX)
      MPI_Allreduce(aux, result, n, MPI_CMPLX, MPI_CMPLX_SUM, comm);
#else
      MPI_Allreduce(aux, result, n, MPI_DOUBLE, MPI_SUM, comm);
#endif    
    PARMS_FREE(aux); 
  }
  return 0;
}
