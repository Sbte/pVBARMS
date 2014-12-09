#include "./include/parms_comm_impl.h"
#include "./include/parms_pc_impl.h"
#include "./include/parms_opt_impl.h"

typedef struct ras_data {
    parms_Operator op;
    parms_Comm handler;
    parms_Map  is;
    FLOAT *rbuf;//changed to buf_recv in the block version
    int nloc;//local size
    int llsize;//point version local size
    int n;//new extended matrix size
    int nn;//new extended matrix size, point version
    int nodv;//number of external variables
    int nsend;//number of sending variables
} *ras_data;

/** Free the memory for struct ras_data.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
static int pc_ras_free(parms_PC *self)
{
    ras_data pc_data;
    parms_Operator op;

    pc_data = (ras_data)(*self)->data;
    op = pc_data->op;
    parms_OperatorFree(&op);
    parms_MapFree(&pc_data->is);
    PARMS_FREE(pc_data);
    (*self)->param->isalloc = false;
    return 0;
}

/** Dump RAS preconditioner.
 *
 *  \param self A preconditioner object.
 *  \param v    A viewer object.
 *  \return 0 on success.
 */
static int pc_ras_view(parms_PC self, parms_Viewer v)
{
    ras_data pc_data;

    pc_data = (ras_data)self->data;
    parms_OperatorView(pc_data->op, v);
    return 0;
}

/** Set up RAS preconditioner.
 *
 *  \param self  A preconditioner object.
 *  \return 0 on success.
 */
static int pc_ras_setup(parms_PC self)
{
    ras_data pc_data;
    parms_Mat A;
    parms_Comm handler;
    parms_FactParam param;
    parms_Operator op;
    void *lmat, *mat_ext;

    /* get the matrix */
    A =  self->A;
    /* get communication handler */
    parms_MatGetCommHandler(A, &handler);

    pc_data = (ras_data)self->data;
    pc_data->nloc = A->is->lsize;
    pc_data->handler = handler;
    pc_data->nodv = parms_CommGetNumRecv(handler);
    pc_data->nsend = parms_CommGetNumSend(handler);
    parms_CommGetRecvBuf(handler, &pc_data->rbuf);

    /* get the local matrix */
    parms_MatGetSubMat(A, &lmat);

    //  parms_vcsr localmat;
    //  localmat = (parms_vcsr)lmat;

    // outputcsmatpa(localmat,"local_data",1);
    //   outputvbmatpa(localmat,"localvb_data",1);


    /* extend the local matrix by including the equations correspond to
     the immediate neighbouring variables */
    parms_MatExtend(A, handler, 0, lmat, &pc_data->n, &mat_ext);//block version

    //  parms_vcsr Amat;
    //  Amat = (parms_vcsr)mat_ext;

    //outputcsmatpa(Amat,"ext_data",1);
    //outputvbmatpa(Amat,"extvb_data",1);
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */

    /* parms_bvcsr Amat; */
    /* Amat = (parms_bvcsr)mat_ext; */

    /* outputvbmatpa(Amat,"b_ext_data",1); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */


    parms_MapCreateFromLocal(&pc_data->is, pc_data->n, 0);

    /* set parameters used for ILU factorization */
    param = self->param;
    param->start = 0;
    printf("param->n = pc_data->n = %d\n",pc_data->n);
    param->n = pc_data->n;
    param->schur_start = param->n;
    if (param->ipar[4] != 0) {
        param->ipar[4] = 0;
    }
    if (param->ipar[5] != 0) {
        param->ipar[5] = 0;
    }

    /* ILU factorization */
    parms_PCILU(self, param, mat_ext, &op);

    /* free the memory for the mat_ext */
    if (self->pcilutype != PCARMS) {//??
        parms_MatFreeSubMat(A, mat_ext);
    }


    pc_data->op = op;
    return 0;
}

/** Apply RAS to vector y.
 *
 *  \f$x = self^{-1}y\f$.
 *
 *  \param self  A RAS preconditioner.
 *  \param y     A right-hand-side vector.
 *  \param x     Solution vector.
 *  \return 0 on success.
 */
static int pc_ras_apply(parms_PC self, FLOAT *y, FLOAT *x)
{//x is the returned solution, y is the input
    ras_data pc_data;
    parms_Operator op;
    parms_Comm handler;
    FLOAT *rbuf;
    FLOAT *x_ext;
    FLOAT *y_ext;
    int nodv, nsend, n, nloc;

    pc_data = (ras_data)self->data;
    op = pc_data->op;
    handler = pc_data->handler;
    rbuf = pc_data->rbuf;
    nsend = pc_data->nsend;
    nodv = pc_data->nodv;
    n = pc_data->n;
    nloc = pc_data->nloc;

    y_ext = (FLOAT *)malloc((n+1)*sizeof(FLOAT));
    x_ext = (FLOAT *)malloc((n+1)*sizeof(FLOAT));

    printf("pc_data->n value is %d.\n", pc_data->n);//%f %p %s %c
    printf("pc_data->nloc value is %d.\n", pc_data->nloc);//%f %p %s %c


    /* exchange interface variables */
    if (nsend) {
        parms_CommDataBegin(handler, y, 0);
    }
    /* copy local variables to the extended vector */
    PARMS_MEMCPY(y_ext, y, pc_data->nloc);

    if (nodv) {
        parms_CommDataEnd(handler);
    }
    /* copy received external interface variables to the extended
     vector */
    if (pc_data->n-pc_data->nloc) {
        PARMS_MEMCPY(&y_ext[nloc], rbuf, pc_data->n-pc_data->nloc);
    }

    /* solve the extended linear system */
    parms_OperatorApply(op, y_ext, x_ext);

    /* ignore external variables -RAS */

    PARMS_MEMCPY(x, x_ext, nloc);
    free(x_ext);
    free(y_ext);

    return 0;
}
/** Get the ratio of the number of nonzero entries of the
 *  preconditioning matrix to that of the original matrix.
 *
 *  \param self   A preconditioner.
 *  \param ratio  A pointer to the ratio.
 *  \return 0 on success.
 */
static int pc_ras_getratio(parms_PC self, double *ratio)
{
    ras_data pc_data;
    parms_Operator op;
    int nnz_mat, nnz_pc;
//    int gnnz_mat, gnnz_pc;

    long int gnnz_mat, gnnz_pc;
    long int nnz_mat_l, nnz_pc_l;

    pc_data = (ras_data)self->data;
    op = pc_data->op;
    parms_OperatorGetNnz(op, &nnz_mat, &nnz_pc);
    nnz_mat_l = (long int)nnz_mat;
    nnz_pc_l = (long int)nnz_pc;

//    MPI_Allreduce(&nnz_mat, &gnnz_mat, 1, MPI_INT, MPI_SUM,
//                  MPI_COMM_WORLD);
//    MPI_Allreduce(&nnz_pc, &gnnz_pc, 1, MPI_INT, MPI_SUM,
//                  MPI_COMM_WORLD);
    MPI_Allreduce(&nnz_mat_l, &gnnz_mat, 1, MPI_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    MPI_Allreduce(&nnz_pc_l, &gnnz_pc, 1, MPI_LONG, MPI_SUM,
                  MPI_COMM_WORLD);
    *ratio = (double)gnnz_pc/(double)gnnz_mat;
    return 0;
}


/** Create a RAS preconditioner.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
int parms_PCCreate_RAS(parms_PC self)
{
    ras_data data;

    PARMS_NEW(data);
    self->data = data;
    self->ops->pc_free     = pc_ras_free;
    self->ops->pc_view     = pc_ras_view;
    self->ops->apply    = pc_ras_apply;
    self->ops->setup    = pc_ras_setup;
    self->ops->getratio = pc_ras_getratio;
    return 0;
}


/** Set up RAS preconditioner.
 *
 *  \param self  A preconditioner object.
 *  \return 0 on success.
 */
static int pc_ras_b_setup(parms_PC self)
{
    ras_data pc_data;
    parms_Mat A;
    parms_Comm handler;
    parms_FactParam param;
    parms_Operator op;
    void *lmat, *mat_ext;

    /* get the matrix */
    A =  self->A;
    /* get communication handler */
    parms_MatGetCommHandler(A, &handler);

    pc_data = (ras_data)self->data;
    pc_data->nloc = A->is->lsize;
    pc_data->llsize = A->is->llsize;//new
    pc_data->handler = handler;
    pc_data->nodv = parms_CommGetNumRecv(handler);
    pc_data->nsend = parms_CommGetNumSend(handler);
    parms_CommGetRecvBuf_b(handler, &pc_data->rbuf);

    /* get the local matrix */
    parms_MatGetSubMat(A, &lmat);

    //  parms_bvcsr localmat;
    //  localmat = (parms_bvcsr)lmat;

    // outputcsmatpa(localmat,"local_data",1);
    //outputvbmatpa(localmat,"localvb_data",1);


    /* extend the local matrix by including the equations correspond to
     the immediate neighbouring variables */
    parms_MatExtend(A, handler, 0, lmat, &pc_data->n, &mat_ext);//block version

    parms_bvcsr Amat;
    Amat = (parms_bvcsr)mat_ext;
    pc_data->nn = Amat->bsz[pc_data->n];
    //outputcsmatpa(Amat,"ext_data",1);
    //outputvbmatpa(Amat,"extvb_data",1);
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */

    /* parms_bvcsr Amat; */
    /* Amat = (parms_bvcsr)mat_ext; */

    /* outputvbmatpa(Amat,"b_ext_data",1); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */


    parms_MapCreateFromLocal(&pc_data->is, pc_data->n, 0);

    /* set parameters used for ILU factorization */
    param = self->param;
    param->start = 0;
    printf("param->n = pc_data->n = %d\n",pc_data->n);
    param->n = pc_data->n;
    param->schur_start = param->n;
    if (param->ipar[4] != 0) {
        param->ipar[4] = 0;
    }
    if (param->ipar[5] != 0) {
        param->ipar[5] = 0;
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    //exit(1);

    /* ILU factorization */
    parms_PCILU(self, param, mat_ext, &op);



    /* free the memory for the mat_ext */
    if (self->pcilutype != PCVBARMS && self->pcilutype != PCVBARMSOLD) {
        parms_MatFreeSubMat(A, mat_ext);
    }


    pc_data->op = op;
    return 0;
}

/** Apply RAS to vector y.
 *
 *  \f$x = self^{-1}y\f$.
 *
 *  \param self  A RAS preconditioner.
 *  \param y     A right-hand-side vector.
 *  \param x     Solution vector.
 *  \return 0 on success.
 */
static int pc_ras_b_apply(parms_PC self, FLOAT *y, FLOAT *x)
{
    ras_data pc_data;
    parms_Operator op;
    parms_Comm handler;
    FLOAT *rbuf;
    FLOAT *x_ext;
    FLOAT *y_ext;
    int nodv, nsend, nn, llsize;//nloc,

    pc_data = (ras_data)self->data;
    op = pc_data->op;
    handler = pc_data->handler;
    rbuf = pc_data->rbuf;
    nsend = pc_data->nsend;
    nodv = pc_data->nodv;
    //n = pc_data->n;//change to point version
    nn = pc_data->nn;//change to point version, extmat size
    //nloc = pc_data->nloc;
    llsize = pc_data->llsize;//here we use llsize//llsize is the point version

    y_ext = (FLOAT *)malloc((nn+1)*sizeof(FLOAT));
    x_ext = (FLOAT *)malloc((nn+1)*sizeof(FLOAT));

    /* exchange interface variables */
    if (nsend) {
        parms_CommDataBegin_b(handler, y, 0);
    }
    /* copy local variables to the extended vector */
    //PARMS_MEMCPY(y_ext, y, pc_data->nloc);//llsize
    PARMS_MEMCPY(y_ext, y, pc_data->llsize);//llsize

    if (nodv) {
        parms_CommDataEnd(handler);
    }
    /* copy received external interface variables to the extended
     vector */
    if (pc_data->n-pc_data->nloc) {
        //PARMS_MEMCPY(&y_ext[nloc], rbuf, pc_data->n-pc_data->nloc);//change the two
        PARMS_MEMCPY(&y_ext[llsize], rbuf, pc_data->nn-pc_data->llsize);//change the two
    }

    /* solve the extended linear system */
    parms_OperatorApply(op, y_ext, x_ext);

    /* ignore external variables -RAS */
    //PARMS_MEMCPY(x, x_ext, nloc);//change nloc
    PARMS_MEMCPY(x, x_ext, llsize);//change nloc
    free(x_ext);
    free(y_ext);

    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */

    return 0;
}


/** Create a RAS preconditioner.
 *
 *  \param self A preconditioner object.
 *  \return 0 on success.
 */
int parms_PCCreate_b_RAS(parms_PC self)
{
    ras_data data;

    PARMS_NEW(data);
    self->data = data;
    self->ops->pc_free     = pc_ras_free;
    self->ops->pc_view     = pc_ras_view;
    self->ops->apply    = pc_ras_b_apply;
    self->ops->setup    = pc_ras_b_setup;
    self->ops->getratio = pc_ras_getratio;
    return 0;
}
