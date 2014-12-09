#include "include/parms_opt_impl.h"
#include "include/parms_comm_impl.h"
#include "DDPQ/protos.h"//new
/*! \typedef parms_dvcsr
  \brief_dvcsr is a synonym for struct parms_dvcsr.
 */
/*! \struct parms_dvcsr
  \brief
 */
typedef struct parms_dvcsr {
    /*! \var diag_mat The diagonal matrix stored in vcsr format.
   */
    parms_vcsr    diag_mat;
    /*!  \var offd_mat The off-diagonal matrix stored in vcsr format.
  */
    parms_vcsr    offd_mat;
    /*! \var mvhandler The parms_Comm object for the matrix-vector product.
   */
    parms_bvcsr   b_diag_mat;//vbsptr
    parms_bvcsr   b_offd_mat;//vbsptr
    parms_Comm    mvhandler;

} *parms_dvcsr;

static int MatVec_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y) 
{
    int lsize, i, j, index, *pj, length;
    parms_dvcsr data;
    parms_Comm  handler;
    parms_vcsr diag_mat, offd_mat;
    parms_Map is;
    FLOAT *pa, *offsetptr, s;

    parms_VecPerm(x, self->is);
    parms_VecPerm(y, self->is);
    //printf("y[0]=%20.16e\n",y[0]);
    /* extract diagonal and off-diagonal matrix */
    data     = (parms_dvcsr)self->data;
    handler  = data->mvhandler;
    diag_mat = data->diag_mat;
    offd_mat = data->offd_mat;
    is       = self->is;
    lsize    = is->lsize;

    parms_CommDataBegin(handler, x, 0);

    /* perform local matrix vector product */
    for (i = 0; i < lsize; i++) {
        s = 0.0;
        length = diag_mat->nnzrow[i];
        pj = diag_mat->pj[i];
        pa = diag_mat->pa[i];
        for (j = 0; j < length; j++) {
            index = pj[j];
            s    += pa[j] * x[index];
        }
        y[i] = s;
    }

    /* off-diagonal part matrix-vector product */
    parms_CommDataEnd(handler);
    offsetptr = handler->buf_recv - lsize;// because in the offd_mat the stating column index is above lsize.


    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        pj     = offd_mat->pj[i];
        pa     = offd_mat->pa[i];
        s      = y[i+is->nint];// the stating row index of offd_mat is 0, so we need to plus the number of internal nodes

        for (j = 0; j < length; j++){
            s   += pa[j] * offsetptr[pj[j]];
            //printf("i = %d, j = %d, pj[j] = %d, length = %d \n", i, j , pj[j], length);getchar();
        }
        y[i+is->nint] = s;
    }
    parms_VecInvPerm(x, self->is);
    parms_VecInvPerm(y, self->is);

    return 0;
}

/** 
 * Setup the matrix
 *
 * This is the most important function in the package. The local
 * variables are divided into two categories: interior variables (not
 * receiving data from other processors) and interface variables
 * (receiving data from other processors). The interior variables are
 * divided into independent variables (not sending data to other
 * processors) and output variables (sending data to other processors)
 * for *asymmetric* pattern. The independent variables are listed first,
 * followd by output variables, and interface variables.  The member
 * variable schur_start is the number of independent
 * variables. i.e., the start position of "Schur" complement in the
 * local linear system.  communication handler for the matrix-vector
 * product is also created.
 *
 * @param self A matrix object.
 *
 * @return 0 on success.
 */
static int MatSetup_dvcsr(parms_Mat self)
{
    parms_dvcsr data;
    parms_vcsr  diag_mat = NULL, offd_mat = NULL;
    parms_vcsr aux_data;
    parms_Map  is;
    parms_Comm handler;
    int i, j, *perm, *iperm, mypid, pid, npro, lsize, k1, k2;
    int *rja=NULL, *tja=NULL, nnz, maxnnz, start, end, nrow;
    int gnodv, cindex, index, gindex;
    int *ginfvars, *ginfv2p, *ginfv2ps, *odperm, *neword, *clp, *info,
            *disp, *pj;
    int npsendrecv, pos, *procslist=NULL;
    int *ptrvsend, *ptrvrecv=NULL, *cmap=NULL;
    FLOAT *ra=NULL, *tma=NULL, val;
    MPI_Comm comm;

    /* extract information arrays */
    is            = self->is;
    npro          = is->npro;
    mypid         = is->pid;
    self->issetup = true;
    data          = (parms_dvcsr)self->data;
    aux_data      = self->aux_data;
    comm          = is->comm;
    handler       = data->mvhandler;
    lsize         = parms_MapGetLocalSize(is);
    perm          = is->perm;
    iperm         = is->iperm;

    //printf("before bsz permutation, lsize = %d\n",lsize);getchar();
    /* free member space */
    if(self->aux_data->space)
        PARMS_FREE(self->aux_data->space);

    /* allocate memory for the diagonal matrix */
    if(self->isreset)
    {
        if(self->resetpattern == SAME_NONZERO_STRUCTURE)
        {
            /* copy data into parallel matrix struct */
            offd_mat = data->offd_mat;
            diag_mat = data->diag_mat;
            /* diagonal part */
            for(i=0; i<diag_mat->n; i++)
            {
                PARMS_MEMCPY(diag_mat->pa[i],  aux_data->pa[i], diag_mat->nnzrow[i]);
            }
            /* off diagonal part */
            for(i=0; i<offd_mat->n; i++)
            {
                index = i+is->nint;
                nrow = diag_mat->nnzrow[index];
                PARMS_MEMCPY(offd_mat->pa[i], &aux_data->pa[index][nrow], offd_mat->nnzrow[i]);
            }
            return 0;
        }
        else
        {
            /* free offd_mat for previous data. No need to free diag_mat
        * for previous data, since array size (lsize) does not change
        */
            diag_mat = data->diag_mat;
            offd_mat = data->offd_mat;
            if(offd_mat->n){
                PARMS_FREE(offd_mat->nnzrow);
                PARMS_FREE(offd_mat->pa);
                PARMS_FREE(offd_mat->pj);
            }
            PARMS_FREE(offd_mat);

            /* free previous communication handler */
            parms_CommFree(&data->mvhandler);
            /* now create new communication handler */
            parms_Comm new_handler;
            parms_CommCreate(&new_handler, comm);
            handler = data->mvhandler = new_handler;
            /* free processor info and table data */
            if(is->vstable) parms_TableFree(&is->vstable);
            if(is->vsend) PARMS_FREE(is->vsend);
        }
    }
    else{
        PARMS_NEW(data->diag_mat);
        diag_mat = data->diag_mat;
        PARMS_NEWARRAY(diag_mat->nnzrow, lsize);
        PARMS_NEWARRAY(diag_mat->pj,     lsize);
        PARMS_NEWARRAY(diag_mat->pa,     lsize);
        diag_mat->n = lsize;
    }

    is->nint = lsize-is->ninf;
    is->schur_start = is->nint;
    k1 = 0;
    k2 = is->nint;
    //printf("before perm");
    for (i = 0; i < lsize; i++) {
        if (perm[i] == -1) {
            perm[i] = k1++;
        }
        else {
            perm[i] = k2++;
        }
    }

    /* set up inverse permutation array */
    for (i = 0; i < lsize; i++) {
        k1 = perm[i];
        iperm[k1] = i;
    }
    is->isperm = true;

    /* permute the rows of the local matrix *///permuted diag_mat first, then assign to aux_data, so permuted both
    for (i = 0; i < lsize; i++) {
        if (aux_data->nnzrow[i]) {
            PARMS_RESIZE(aux_data->pj[i], aux_data->nnzrow[i]);
            PARMS_RESIZE(aux_data->pa[i], aux_data->nnzrow[i]);
        }
        k1 = perm[i];
        diag_mat->pj[k1]     = aux_data->pj[i];
        diag_mat->pa[k1]     = aux_data->pa[i];
        diag_mat->nnzrow[k1] = aux_data->nnzrow[i];
    }

    PARMS_MEMCPY(aux_data->pj,     diag_mat->pj,     lsize);
    PARMS_MEMCPY(aux_data->pa,     diag_mat->pa,     lsize);
    PARMS_MEMCPY(aux_data->nnzrow, diag_mat->nnzrow, lsize);

    /* permute the column indices of the matrix which are corresponding
     to the interior variables *///only aux_data
    for (i = 0; i < is->nint; i++) {
        rja = aux_data->pj[i];
        nnz = aux_data->nnzrow[i];
        ra  = aux_data->pa[i];
        for (j = 0; j < nnz; j++) {
            rja[j] = perm[rja[j]];	/* in C style */
        }
    }

    /* permute the local column indices of the matrix which are
     corresponding to the interface variables *///also aux_data
    maxnnz = 0;
    for (i = is->nint; i < lsize; i++) {
        rja = aux_data->pj[i];
        nnz = aux_data->nnzrow[i];
        if (maxnnz < nnz) {
            maxnnz =  nnz;
        }
        diag_mat->nnzrow[i] = 0;
        for (j = 0; j < nnz; j++) {
            if (rja[j] >= 0) {
                diag_mat->nnzrow[i]++;
                rja[j] = perm[rja[j]]; /* in C style */
            }
        }
    }

    if (maxnnz) {
        PARMS_NEWARRAY(tja, maxnnz);
        PARMS_NEWARRAY(tma, maxnnz);
    }

    /* allocate memory for the off-diagonal matrix */
    PARMS_NEW(data->offd_mat);
    offd_mat = data->offd_mat;
    offd_mat->n = is->ninf;
    if (is->ninf) {
        PARMS_NEWARRAY(offd_mat->nnzrow, is->ninf);
        PARMS_NEWARRAY(offd_mat->pj,     is->ninf);
        PARMS_NEWARRAY(offd_mat->pa,     is->ninf);
    }
    handler->nodv = parms_TableGetSize(is->table) - lsize;
    if (handler->nodv) {
        PARMS_NEWARRAY(handler->odvlist, handler->nodv);
    }

    /* list column indices of diagonal part first */
    for (i = is->nint; i < lsize; i++) {
        rja    = aux_data->pj[i];
        ra     = aux_data->pa[i];
        nnz    = aux_data->nnzrow[i];
        index  = i-is->nint;
        k1     = 0;
        k2     = diag_mat->nnzrow[i];
        nrow   = diag_mat->nnzrow[i];
        for (j = 0; j < nnz; j++) {
            if (rja[j] >= 0) {
                tja[k1]   = rja[j];
                tma[k1++] = ra[j];
            }
            else {
                tja[k2]   = -rja[j]-1;
                clp = parms_MapGlobalToLocal(is, tja[k2]);
                if (clp != NULL) {
                    handler->odvlist[*clp-lsize] = tja[k2];
                }
                tma[k2++] = ra[j];
            }
        }
        PARMS_MEMCPY(rja, tja, nnz);//(dse, sou, nnz)
        PARMS_MEMCPY(ra,  tma, nnz);
        offd_mat->pj[index]     = &rja[nrow];
        offd_mat->pa[index]     = &ra[nrow];
        offd_mat->nnzrow[index] = nnz-nrow;
    }
    if (maxnnz) {
        PARMS_FREE(tja);
        PARMS_FREE(tma);
    }

    /* sort the column indices of the diagonal matrix in ascending
     order */
    for (i = 0; i < lsize; i++) {
        nnz = diag_mat->nnzrow[i];
        rja = diag_mat->pj[i];
        ra  = diag_mat->pa[i];
        for (j = 0; j < nnz-1; j++) {
            for (k1 = j+1; k1 < nnz; k1++) {
                if (rja[j] > rja[k1]) {
                    index = rja[j];
                    rja[j] = rja[k1];
                    rja[k1] = index;

                    val = ra[j];
                    ra[j] = ra[k1];
                    ra[k1] = val;
                }
            }
        }
    }

    PARMS_NEWARRAY(info, npro);
    PARMS_NEWARRAY(disp, npro+1);

    MPI_Allgather(&handler->nodv, 1, MPI_INT, info, 1, MPI_INT, comm);
    disp[0] = 0;
    for (i = 0; i < npro; i++) {
        disp[i+1] = disp[i] + info[i];
    }

    /* total number of interface variables */
    gnodv = disp[npro];
    if (gnodv) {
        PARMS_NEWARRAY(ginfvars,  gnodv);
        PARMS_NEWARRAY0(ginfv2p,  gnodv);
        PARMS_NEWARRAY(ginfv2ps,  gnodv);
        MPI_Allgatherv(handler->odvlist, handler->nodv, MPI_INT, ginfvars,
                       info, disp, MPI_INT, comm);

        /* variables sent to other processors */
        for (i = 0; i < npro; i++) {
            if (i != mypid) {
                for (j = disp[i]; j < disp[i+1]; j++) {
                    clp = parms_MapGlobalToLocal(is, ginfvars[j]);
                    if (clp != NULL && *clp < lsize) {
                        ginfv2p[j] = mypid;
                    }
                }
            }
        }

        MPI_Allreduce(ginfv2p, ginfv2ps, gnodv, MPI_INT, MPI_SUM, comm);
        PARMS_FREE(ginfv2p);
        PARMS_NEWARRAY(handler->procs_recv, npro);
        for (i = 0; i < npro; i++) {
            info[i] = 0;
        }
        /* processor IDs from which processor mypid receives data */
        handler->nprecv = 0;
        for (i = disp[mypid]; i < disp[mypid+1]; i++) {
            pid = ginfv2ps[i];
            if (info[pid] == 0) {
                handler->procs_recv[handler->nprecv++] = pid;
            }
            info[pid]++;
        }

        if (handler->nprecv) {
            PARMS_RESIZE(handler->procs_recv, handler->nprecv);

            /* sort neighbouring processors received data from in ascending order */
            for (i = 0; i < handler->nprecv-1; i++) {
                for (j = i+1; j < handler->nprecv; j++) {
                    if (handler->procs_recv[i] > handler->procs_recv[j]) {
                        pid = handler->procs_recv[i];
                        handler->procs_recv[i] = handler->procs_recv[j];
                        handler->procs_recv[j] = pid;
                    }
                }
            }
        }
        else {
            PARMS_FREE(handler->procs_recv);
        }

        PARMS_NEWARRAY(handler->ptrvrecv, handler->nprecv+1);
        handler->ptrvrecv[0] = 0;
        for (i = 0; i < handler->nprecv; i++) {
            pid = handler->procs_recv[i];
            handler->ptrvrecv[i+1] = handler->ptrvrecv[i] + info[pid];
            info[pid] = i;
        }

        if (handler->nodv) {
            PARMS_NEWARRAY(odperm, handler->nodv);
            PARMS_NEWARRAY(handler->buf_recv, handler->nodv);

            /* reorder external variables: the data received from the same
     processor are listed consecutively, processor by processor */
            j = disp[mypid];
            for (i = 0; i < handler->nodv; i++) {
                pid = ginfv2ps[i+j];
                index = info[pid];
                odperm[i] = handler->ptrvrecv[index]++;
            }
            for (i = handler->nprecv-1; i > 0; i--) {
                handler->ptrvrecv[i] = handler->ptrvrecv[i-1];
            }
            handler->ptrvrecv[0] = 0;

            /* reorder the external interface variables */
            PARMS_NEWARRAY(neword, handler->nodv);
            for (i = 0; i < handler->nodv; i++) {
                index = odperm[i];
                gindex = handler->odvlist[i];
                parms_TablePut(is->table, gindex, index+lsize);// insert the pair (key, value) into the table.
                neword[index] = gindex;
            }

            PARMS_FREE(handler->odvlist);

            PARMS_FREE(odperm);
            handler->odvlist = neword;
        }
        /* change the global column indices of entries in external
       matrix to the local ones */
        //outputcsmatpa(offd_mat,"offd_datainsetupbefore_local",1);

        for (i = 0; i < offd_mat->n; i++) {
            nnz = offd_mat->nnzrow[i];
            pj  = offd_mat->pj[i];
            for (j = 0; j < nnz; j++) {
                cindex = pj[j];
                clp  = parms_MapGlobalToLocal(is, cindex);
                pj[j]  = *clp;

                //if (pid == 0 ) printf("row = %d, cindex = %d, clp = %d \n",i, cindex,*clp);//getchar();

            }
        }
        PARMS_FREE(ginfv2ps);
        PARMS_FREE(ginfvars);
    }
    else {
        handler->nprecv = 0;
        PARMS_NEWARRAY(handler->ptrvrecv, 1);
        handler->ptrvrecv[0] = 0;
    }

    /* set up send information */
    MPI_Allgather(&handler->nprecv, 1, MPI_INT, info, 1, MPI_INT, comm);
    disp[0] = 0;
    for (i = 0; i < npro; i++) {
        disp[i+1] = disp[i] + info[i];
    }

    if (disp[npro]) {
        PARMS_NEWARRAY(procslist, disp[npro]);
        MPI_Allgatherv(handler->procs_recv, handler->nprecv, MPI_INT,
                       procslist, info, disp, MPI_INT, comm);
    }

    handler->npsend = 0;
    PARMS_NEWARRAY(handler->procs_send, npro);
    for (i = 0; i < npro; i++) {
        if (i != mypid) {
            for (j = disp[i]; j < disp[i+1]; j++) {
                if (procslist[j] == mypid) {
                    handler->procs_send[handler->npsend++] = i;
                    break;
                }
            }
        }
    }

    if (disp[npro]) {
        PARMS_FREE(procslist);
    }

    if (handler->npsend) {
        PARMS_RESIZE(handler->procs_send, handler->npsend);
    }
    else {
        PARMS_FREE(handler->procs_send);
    }

    PARMS_NEWARRAY(handler->ptrvsend, handler->npsend+1);
    handler->ptrvsend[0] = 0;
    ptrvsend = handler->ptrvsend;
    ptrvrecv = handler->ptrvrecv;
    npsendrecv = handler->npsend + handler->nprecv;
    if (npsendrecv) {
        PARMS_NEWARRAY(handler->status_send, npsendrecv);
        PARMS_NEWARRAY(handler->req_send,    npsendrecv);
        handler->status_recv = handler->status_send + handler->npsend;
        handler->req_recv    = handler->req_send + handler->npsend;
        for (i = 0; i < handler->npsend; i++) {
            MPI_Irecv(&disp[i], 1, MPI_INT, handler->procs_send[i],
                      100, comm, &handler->req_send[i]);
        }

        for (i = 0; i < handler->nprecv; i++) {
            info[i] = ptrvrecv[i+1] - ptrvrecv[i];
            MPI_Isend(&info[i], 1, MPI_INT, handler->procs_recv[i],
                      100, comm, &handler->req_recv[i]);
        }

        MPI_Waitall(handler->nprecv, handler->req_recv,
                    handler->status_recv);
        MPI_Waitall(handler->npsend, handler->req_send,
                    handler->status_send);

        for (i = 0; i < handler->npsend; i++) {
            ptrvsend[i+1] = ptrvsend[i] + disp[i];
        }

        nnz = ptrvsend[handler->npsend];
        if (nnz) {
            PARMS_NEWARRAY(handler->buf_send,   nnz);
            PARMS_NEWARRAY(handler->vlist_send, nnz);
        }

        handler->mdata_send = 0;
        if (handler->npsend) {
            for (i = 0; i < handler->npsend; i++) {
                nnz = ptrvsend[i+1] - ptrvsend[i];
                if (handler->mdata_send < nnz) {
                    handler->mdata_send = nnz;
                }
            }
        }

        /* handler->vlist_send stores the variables to be sent to adjacent
       processors in global labeling */
        for (i = 0; i < handler->npsend; i++) {
            pos = ptrvsend[i];
            nnz = ptrvsend[i+1] - ptrvsend[i];
            MPI_Irecv(&handler->vlist_send[pos], nnz, MPI_INT,
                      handler->procs_send[i], 200, comm,
                      &handler->req_send[i]);
        }

        for (i = 0; i < handler->nprecv; i++) {
            pos = ptrvrecv[i];
            nnz = ptrvrecv[i+1] - ptrvrecv[i];
            MPI_Isend(&handler->odvlist[pos], nnz, MPI_INT,
                      handler->procs_recv[i], 200, comm,
                      &handler->req_recv[i]);
        }

        MPI_Waitall(handler->nprecv, handler->req_recv,
                    handler->status_recv);
        MPI_Waitall(handler->npsend, handler->req_send,
                    handler->status_send);

        /* change global labelings  to local ones */
        is->ninf_send = 0;
        parms_TableCreate(&is->vstable, NULL, is->ninf);
        if (is->isperm) {
            for (i = 0; i < handler->npsend; i++) {
                for (j = ptrvsend[i]; j < ptrvsend[i+1]; j++) {
                    cindex = handler->vlist_send[j];
                    clp  = parms_MapGlobalToLocal(is, cindex);
                    index  = perm[*clp];
                    handler->vlist_send[j] = index;
                    clp = parms_TableGet(is->vstable, index);
                    if (clp == NULL) {
                        parms_TablePut(is->vstable, index, is->ninf_send++);
                    }
                }
            }
        }

        /* vsend is an array to store variables sent to other processors */
        if (is->ninf_send) {
            PARMS_NEWARRAY(is->vsend, is->ninf_send);
            if (is->nint) {
                PARMS_NEWARRAY(cmap, is->nint);
            }
            for (i = 0; i < is->nint; i++) {
                cmap[i] = -1;
            }

            for (i = 0; i < handler->npsend; i++) {
                for (j = ptrvsend[i]; j < ptrvsend[i+1]; j++) {
                    index = handler->vlist_send[j];
                    clp = parms_TableGet(is->vstable, index);
                    is->vsend[*clp] = index;

                }
            }
            for (i = 0; i < is->ninf_send; i++) {
                if (is->vsend[i] < is->nint) {
                    cmap[is->vsend[i]] = --is->schur_start;
                }
            }

            /* reorder the local variables in the following order:
       * 1. independent variables which are uncoupled with variables in
       *    adjacent processors.
       * 2. variables sent to adjacent processors but the
       *    corresponding rows are decoupled with variables in
       *    the adjacent processors. (asymmetric pattern)
       * 3. interface variables - coupled with variables in adjacent
       *    processors in the sense of receiving data from
       *    the adjacent processors.
       */
            if (is->schur_start != is->nint) {
                k1 = 0;
                for (i = 0; i < is->nint; i++) {
                    if (cmap[i] == -1) {
                        cmap[i] = k1++;
                    }
                }
                /* permute column indices of the local matrix */
                for (i = 0; i < lsize; i++) {
                    nnz = diag_mat->nnzrow[i];
                    pj = diag_mat->pj[i];
                    for (j = 0; j < nnz; j++) {
                        if (pj[j] < is->nint) {
                            pj[j] = cmap[pj[j]];
                        }
                    }
                }

                /* permute rows which are corresponding to interior variables */
                for (i = 0; i < is->nint; i++) {
                    k1 = cmap[i];
                    diag_mat->nnzrow[k1] = aux_data->nnzrow[i];
                    diag_mat->pj[k1]     = aux_data->pj[i];
                    diag_mat->pa[k1]     = aux_data->pa[i];
                }
                PARMS_MEMCPY(aux_data->nnzrow, diag_mat->nnzrow, is->nint);
                PARMS_MEMCPY(aux_data->pj,     diag_mat->pj,     is->nint);
                PARMS_MEMCPY(aux_data->pa,     diag_mat->pa,     is->nint);

                /* update perm and iperm */
                for (i = 0; i < lsize; i++) {
                    if (perm[i] < is->nint) {
                        perm[i] = cmap[perm[i]];
                    }
                }

                for (i = 0; i < lsize; i++) {
                    k1 = perm[i];
                    iperm[k1] = i;
                }

                /* update vlist_send */
                for (i = 0; i < handler->npsend; i++) {
                    start = handler->ptrvsend[i];
                    end   = handler->ptrvsend[i+1];
                    for (j = start; j < end; j++) {
                        if (handler->vlist_send[j] < is->nint) {
                            handler->vlist_send[j] = cmap[handler->vlist_send[j]];
                        }
                    }
                }

                /* update vsend */
                for (i = 0; i < is->ninf_send; i++) {
                    if (is->vsend[i] < is->nint) {
                        is->vsend[i] = cmap[is->vsend[i]];
                        parms_TablePut(is->vstable, is->vsend[i], i);
                    }
                }
            }
            if (is->nint) {
                PARMS_FREE(cmap);
            }

        }

    }
    /* outputcsmatpa(diag_mat,"diag_data",1); */
    /* outputcsmatpa(offd_mat,"offd_data",1); */

    //output_intvectorpa("columnglobal",handler->odvlist,0, handler->nodv);
    PARMS_FREE(info);
    PARMS_FREE(disp);

    return 0;
}  


static int MatSetCommType_dvcsr(parms_Mat self, COMMTYPE ctype)
{
    parms_dvcsr data;
    parms_Comm  handler;

    data           = (parms_dvcsr)self->data;
    handler        = data->mvhandler;
    handler->ctype = ctype;
    return 0;
}

static int MatMVPY_dvcsr(parms_Mat self, FLOAT alpha, FLOAT *x,
                         FLOAT beta, FLOAT *y, FLOAT *z)
{
    int lsize, i, j, index, *pj, length;
    parms_dvcsr data;
    parms_Comm handler;
    parms_vcsr diag_mat, offd_mat;
    parms_Map is;
    FLOAT *pa, *offsetptr, s;

    /* no need to permute z here, since it only stores the result */
    parms_VecPerm(x, self->is);
    parms_VecPerm(y, self->is);

    /* extract diagonal and off-diagonal matrix */
    data     = (parms_dvcsr)self->data;
    handler  = data->mvhandler;
    diag_mat = data->diag_mat;
    offd_mat = data->offd_mat;
    is       = self->is;
    lsize    = is->lsize;

    /* post receive and send actions */
    parms_CommDataBegin(handler, x, 0);

    /* perform local matrix vector product */
    for (i = 0; i < lsize; i++) {
        s = beta * y[i];
        length = diag_mat->nnzrow[i];
        pj = diag_mat->pj[i];
        pa = diag_mat->pa[i];
        for (j = 0; j < length; j++) {
            index = pj[j];
            s    += alpha * pa[j] * x[index];
        }
        z[i] = s;
    }

    parms_CommDataEnd(handler);

    offsetptr = handler->buf_recv - lsize;
    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        pj     = offd_mat->pj[i];
        pa     = offd_mat->pa[i];
        s      = z[i+is->nint];
        for (j = 0; j < length; j++) {
            s   += alpha * pa[j] * offsetptr[pj[j]];
        }
        z[i+is->nint] = s;
    }

    /* inverse permutation. z now has permutation of y, so has to be inverted as well */
    parms_VecInvPerm(x, self->is);
    parms_VecInvPerm(y, self->is);
    parms_VecInvPerm(z, self->is);
    return 0;
}

static int MatVecSchur_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y, int
                             pos)
{
    parms_dvcsr data;
    parms_Comm  handler;
    parms_Map   is;
    parms_vcsr  offd_mat;
    FLOAT       *offsetptr, *pa, s;
    int i, j, length, *pj, lsize;

    is       = self->is;
    lsize    = is->lsize;
    data     = (parms_dvcsr)self->data;
    handler  = data->mvhandler;
    offd_mat = data->offd_mat;

    if ((x == NULL) || (y == NULL)) {
        return 0;
    }

    parms_CommDataBegin(handler, x, pos);
    for (i = pos; i < is->nint; i++) {
        y[i-pos] = 0.0;
    }
    parms_CommDataEnd(handler);
    offsetptr = handler->buf_recv - lsize;
    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        pj     = offd_mat->pj[i];
        pa     = offd_mat->pa[i];
        s      = 0.0;
        for (j = 0; j < length; j++) {
            s   += pa[j] * offsetptr[pj[j]];
        }
        y[i+is->nint-pos] = s;
    }
    return 0;
}

static int MatGetDiag_dvcsr(parms_Mat self, void **mat)
{
    parms_dvcsr data;
    parms_vcsr  diag_mat;

    data = (parms_dvcsr)self->data;
    diag_mat = data->diag_mat;
    //printf("Before if judge");
    if (self->ilutype == PCARMS) {
        parms_vcsr diag;
        int i, nnz;
        //printf("in MatGetDiag_dvcsr a new diag");
        PARMS_NEW(diag);
        diag->n = diag_mat->n;
        PARMS_NEWARRAY(diag->nnzrow, diag->n);
        PARMS_NEWARRAY(diag->pj,     diag->n);
        PARMS_NEWARRAY(diag->pa,     diag->n);
        for (i = 0; i < diag->n; i++) {
            nnz = diag->nnzrow[i] = diag_mat->nnzrow[i];
            if (nnz) {
                PARMS_NEWARRAY(diag->pj[i], nnz);
                PARMS_NEWARRAY(diag->pa[i], nnz);
                PARMS_MEMCPY(diag->pj[i], diag_mat->pj[i], nnz);
                PARMS_MEMCPY(diag->pa[i], diag_mat->pa[i], nnz);
            }
        }
        *mat = diag;
        return 0;
    }
    *mat = diag_mat;
    return 0;
}

static int MatGetSubMat_dvcsr(parms_Mat self, void **mat)
{
    *mat = self->aux_data;
    return 0;
}

static int MatExtend_dvcsr(parms_Mat self, parms_Comm handler, int
                           start, void *mat, int *nn, void
                           **ext_mat)
{
    parms_vcsr lmat, new_mat;
    parms_Map  is;
    MPI_Request *sreq=NULL, *rreq=NULL;
    MPI_Status *sstatus=NULL, *rstatus=NULL;
    MPI_Comm comm;
    FLOAT **a_send=NULL, *pa_rbuf=NULL, *pa_sbuf=NULL, *rowa=NULL;
    int **ja_send=NULL;
    int *ia_send=NULL, *pia_send=NULL, *pia_recv=NULL, *pj_rbuf=NULL, *pj_sbuf=NULL;
    int *pindex, *rowj, *rbuf=NULL, *sbuf=NULL;
    int lsize, jcol, i, j, index, nnz, nodv, numsend, numrecv;
    int n, pos, lindex, tag, length, k, m;

//    int reveiverlength = 0;



    lmat = (parms_vcsr)mat;
    PARMS_NEW(new_mat);
    is = self->is;
    lsize = is->lsize;
    outputcsmatpa(lmat,"lmat",1);

    MPI_Comm_dup(handler->comm, &comm);

    //printf("is->ninf_send = %d\n", is->ninf_send);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    printf("is->nint = %d, pid = %d \n",is->nint, myid);
    printf("is->ninf = %d, pid = %d \n",is->ninf, myid);
    printf("is->ninf_send = %d, pid = %d \n",is->ninf_send, myid);
    //printf("The number of stored entries is %d in local block form matrix, pid = %d \n",memVBMat(Amat), myid);


    if (is->ninf_send) {
        PARMS_NEWARRAY(a_send,  is->ninf_send);//2D pointer to values
        PARMS_NEWARRAY(ja_send, is->ninf_send);//2D pointer to index
        PARMS_NEWARRAY(ia_send, is->ninf_send);//array to store nnz

        /* copy the matrix corresponding to the interface variables to
       working arrays. the column indices are stored in global numbering */
        for (i = 0; i < is->ninf_send; i++) {
            index = is->vsend[i] - start;//the index of variables to send out
            nnz = ia_send[i] = lmat->nnzrow[index];
            PARMS_NEWARRAY(a_send[i],  nnz);
            PARMS_NEWARRAY(ja_send[i], nnz);

            //      output_intvectorpa("ia_send",ia_send,0, is->ninf_send);


            nnz = lmat->nnzrow[index];
            for (j = 0; j < nnz; j++) {
                a_send[i][j] = lmat->pa[index][j];
                jcol = lmat->pj[index][j];
                if (jcol < lsize) {
                    if (is->isperm) {
                        jcol = is->iperm[jcol];
                        ja_send[i][j] = is->lvars[jcol];
                    }
                    else {
                        ja_send[i][j] = is->lvars[jcol];
                    }
                }
                else {
                    jcol = handler->odvlist[jcol-lsize];
                    ja_send[i][j] = jcol;
                }
            }
        }
    }


    numsend = handler->ptrvsend[handler->npsend];
    numrecv = handler->ptrvrecv[handler->nprecv];

    nodv = handler->nodv;
    n = lmat->n + nodv;

    new_mat->n = n;
    PARMS_NEWARRAY(new_mat->nnzrow, n);
    PARMS_NEWARRAY(new_mat->pj,     n);
    PARMS_NEWARRAY(new_mat->pa,     n);

    /* copy the local matrix to new_mat */
    for (i = 0; i < lmat->n; i++) {
        nnz = new_mat->nnzrow[i] = lmat->nnzrow[i];
        if (nnz) {
            PARMS_NEWARRAY(new_mat->pj[i], nnz);
            PARMS_NEWARRAY(new_mat->pa[i], nnz);
            PARMS_MEMCPY(new_mat->pj[i], lmat->pj[i], nnz);
            PARMS_MEMCPY(new_mat->pa[i], lmat->pa[i], nnz);
        }
    }

    if (handler->npsend+handler->nprecv) {
        PARMS_NEWARRAY(sstatus, handler->npsend+handler->nprecv);
        PARMS_NEWARRAY(sreq,    handler->npsend+handler->nprecv);
        rstatus = sstatus + handler->npsend;
        rreq    = sreq + handler->npsend;
    }

    tag = 800;

    /* exchange external matrix */
    if (numrecv) {
        PARMS_NEWARRAY(rbuf, numrecv);//rbuf is to store the nnz of the rows
    }

    /* receive the number of entries for each row in the extended
     matrix */
    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        MPI_Irecv(&rbuf[pos], length, MPI_INT, handler->procs_recv[i],
                  tag, comm, &rreq[i]);
    }

    if (numsend) {
        PARMS_NEWARRAY(sbuf, numsend);//sbuf is to store the nnz of the rows
    }

    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);

            /* if(myid==2) */
            /* 	printf("pos+j = %d, lindex = %d, newlindex = %d, pid = %d \n",pos+j, lindex, *pindex ,myid); */

            lindex = *pindex;
            sbuf[pos+j] = ia_send[lindex];
        }
        MPI_Isend(&sbuf[pos], length, MPI_INT, handler->procs_send[i],
                  tag, comm, &sreq[i]);
    }
    MPI_Waitall(handler->nprecv, rreq, rstatus);
    MPI_Waitall(handler->npsend, sreq, sstatus);

    // output_intvectorpa("rbuf",rbuf,0, numrecv);


    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        for (j = 0; j < length; j++) {
            new_mat->nnzrow[lmat->n+pos+j] = rbuf[pos+j];
        }
    }

    j = 0;
    for (i = 0; i < nodv; i++) {
        length = new_mat->nnzrow[i+lsize];
        j += length;
        if (length) {
            PARMS_NEWARRAY(new_mat->pj[i+lmat->n], length);
            PARMS_NEWARRAY(new_mat->pa[i+lmat->n], length);
        }
    }
    if (j) {// j is the sum of the receiver's length
//        reveiverlength = j;
        PARMS_NEWARRAY(pj_rbuf, j);//int index
        PARMS_NEWARRAY(pa_rbuf, j);//FLOAT value, a array for all the rows
    }

    PARMS_NEWARRAY(pia_send, handler->npsend+1);
    PARMS_NEWARRAY(pia_recv, handler->nprecv+1);
    pia_send[0] = 0;
    pia_recv[0] = 0;

    /* receive column indices of each row */
    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            m = new_mat->nnzrow[lmat->n+pos+j];
            k += m;
        }
        pia_recv[i+1] = pia_recv[i] + k;
        MPI_Irecv(&pj_rbuf[pia_recv[i]],  k, MPI_INT,
                handler->procs_recv[i], tag, comm, &rreq[i]);
    }

    k = 0;
    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            k += m;
        }
    }

    if (k) {
        PARMS_NEWARRAY(pj_sbuf, k);
        PARMS_NEWARRAY(pa_sbuf, k);
    }

    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            PARMS_MEMCPY(&pj_sbuf[pia_send[i]+k], ja_send[lindex], m);
            k += m;
        }
        pia_send[i+1] = pia_send[i] + k;
        MPI_Isend(&pj_sbuf[pia_send[i]], k, MPI_INT, handler->procs_send[i],
                tag, comm, &sreq[i]);
    }
    MPI_Waitall(handler->nprecv, rreq, rstatus);
    MPI_Waitall(handler->npsend, sreq, sstatus);

    //output_intvectorpa("pj_rbuf",pj_rbuf,0, reveiverlength);


    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            m = new_mat->nnzrow[lmat->n+pos+j];
            PARMS_MEMCPY(new_mat->pj[lmat->n+pos+j], &pj_rbuf[pia_recv[i]+k],
                    m);
            k += m;
        }
    }

    /* receive values of each row */

#if defined(DBL_CMPLX)
    for (i = 0; i < handler->nprecv; i++) {
        k = pia_recv[i+1] - pia_recv[i];
        MPI_Irecv(&pa_rbuf[pia_recv[i]], k, MPI_CMPLX,
                handler->procs_recv[i], tag, comm, &rreq[i]);
    }

    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            PARMS_MEMCPY(&pa_sbuf[pia_send[i]+k], a_send[lindex], m);
            k += m;
        }
        MPI_Isend(&pa_sbuf[pia_send[i]], k, MPI_CMPLX,
                handler->procs_send[i], tag, comm, &sreq[i]);
    }
#else
    for (i = 0; i < handler->nprecv; i++) {
        k = pia_recv[i+1] - pia_recv[i];
        MPI_Irecv(&pa_rbuf[pia_recv[i]], k, MPI_DOUBLE,
                handler->procs_recv[i], tag, comm, &rreq[i]);
    }

    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            PARMS_MEMCPY(&pa_sbuf[pia_send[i]+k], a_send[lindex], m);
            k += m;
        }
        MPI_Isend(&pa_sbuf[pia_send[i]], k, MPI_DOUBLE,
                handler->procs_send[i], tag, comm, &sreq[i]);
    }
#endif

    MPI_Waitall(handler->nprecv, rreq, rstatus);
    MPI_Waitall(handler->npsend, sreq, sstatus);

    MPI_Comm_free(&comm);
    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            m = new_mat->nnzrow[lmat->n+pos+j];
            PARMS_MEMCPY(new_mat->pa[lmat->n+pos+j], &pa_rbuf[pia_recv[i]+k],
                    m);
            k += m;
        }
    }

    /* free temporary arrays */
    if (handler->npsend) {
        PARMS_FREE(pj_sbuf);
        PARMS_FREE(pa_sbuf);
    }

    if (handler->nprecv) {
        PARMS_FREE(pj_rbuf);
        PARMS_FREE(pa_rbuf);
    }

    if (numrecv) {
        PARMS_FREE(rbuf);
    }

    if (numsend) {
        PARMS_FREE(sbuf);
    }


    PARMS_FREE(pia_recv);
    PARMS_FREE(pia_send);

    if (handler->npsend+handler->nprecv) {
        PARMS_FREE(sreq);
        PARMS_FREE(sstatus);
    }

    if (is->ninf_send) {
        for (i = 0; i < is->ninf_send; i++) {
            if (ia_send[i]) {
                PARMS_FREE(ja_send[i]);
                PARMS_FREE(a_send[i]);
            }
        }

        PARMS_FREE(ia_send);
        PARMS_FREE(ja_send);
        PARMS_FREE(a_send);
    }

    /* change the global column indices into local ones */
    for (i = lmat->n; i < n; i++) {
        nnz = new_mat->nnzrow[i];
        rowj = new_mat->pj[i];
        rowa = new_mat->pa[i];
        k = 0;
        for (j = 0; j < nnz; j++) {
            jcol = rowj[j];
            pindex = parms_TableGet(is->table, jcol);//global to local
            if (pindex) {
                if (*pindex < lsize) {
                    rowj[j] = is->perm[*pindex] - start;
                }
                else if (*pindex >= lsize) {
                    rowj[j] = *pindex - start;
                }
                if (k < j) {
                    rowj[k] = rowj[j];
                    rowa[k] = rowa[j];
                }
                k++;
            }
        }
        new_mat->nnzrow[i] = k;
        if (k) {
            PARMS_RESIZE(new_mat->pj[i], k);
            PARMS_RESIZE(new_mat->pa[i], k);
        }
    }
    *ext_mat = new_mat;
    *nn = n;
    return 0;
}

static int MatFreeSubMat_dvcsr(parms_Mat self, void *mat)
{
    parms_vcsr vmat;
    int i, nnz, isalloc;

    isalloc = self->isalloc;

    if(isalloc){
        vmat = (parms_vcsr)mat;
        for (i = 0; i < vmat->n; i++) {
            nnz = vmat->nnzrow[i];
            if (nnz) {
                PARMS_FREE(vmat->pj[i]);
                PARMS_FREE(vmat->pa[i]);
            }
        }

        if (vmat->n) {
            PARMS_FREE(vmat->nnzrow);
            PARMS_FREE(vmat->pj);
            PARMS_FREE(vmat->pa);
        }

        PARMS_FREE(vmat);

        mat = NULL;
    }

    return 0;
}

static int MatGetHandler_dvcsr(parms_Mat self, parms_Comm *handler)
{
    parms_dvcsr data;

    data = (parms_dvcsr)self->data;
    *handler = data->mvhandler;
    return 0;
}

static struct parms_Mat_ops parms_matops_dvcsr = {//make the global variable static
    MatVec_dvcsr,
    MatSetup_dvcsr,
    MatSetCommType_dvcsr,
    MatMVPY_dvcsr,
    MatGetDiag_dvcsr,
    MatGetSubMat_dvcsr,
    MatExtend_dvcsr,
    MatVecSchur_dvcsr,
    MatFreeSubMat_dvcsr,
    MatGetHandler_dvcsr
};

int parms_MatFree_dvcsr(parms_Mat *self)
{
    parms_dvcsr   data;
    parms_vcsr    diag_mat, offd_mat;

    data = (parms_dvcsr)(*self)->data;

    parms_CommFree(&data->mvhandler);
    diag_mat = data->diag_mat;
    PARMS_FREE(diag_mat->nnzrow);
    PARMS_FREE(diag_mat->pa);
    PARMS_FREE(diag_mat->pj);
    PARMS_FREE(diag_mat);

    offd_mat = data->offd_mat;
    if (offd_mat->n) {
        PARMS_FREE(offd_mat->nnzrow);
        PARMS_FREE(offd_mat->pa);
        PARMS_FREE(offd_mat->pj);
    }
    PARMS_FREE(offd_mat);

    PARMS_FREE(data);
    return 0;
}

int parms_MatView_dvcsr(parms_Mat self, parms_Viewer v)
{
    int i, j, lsize, pid, npro, length;
    int *rja, *lvlist, *iperm, index;
    FLOAT *ra;
    parms_dvcsr data;
    parms_vcsr diag_mat, offd_mat;
    parms_Map is;
    FILE *fp;

    is      = self->is;
    lvlist  = is->lvars;
    iperm   = is->iperm;
    lsize   = is->lsize;
    pid     = is->pid;
    npro    = is->npro;
    data    = (parms_dvcsr)self->data;

    parms_ViewerGetFP(v, &fp);

    if (pid == 0) {
        fprintf(fp, "There are %d processors available\n", npro);
        fprintf(fp, "The format of output local equations is:\n");
        fprintf(fp, "(pid,local_row_index)=>(pid,global_row_index)\n");
        fprintf(fp, "(pid,local_row_index, local_column_index) = value\n");
    }
    fprintf(fp, "\n");

    for (i = 0; i < lsize; i++) {
        index = iperm[i];
        fprintf(fp, "(%d,%d)=>(%d,%d)\n", pid, i, pid, lvlist[index]);
    }
    fprintf(fp, "\n");

    diag_mat = data->diag_mat;
    offd_mat = data->offd_mat;
    fprintf(fp, "Local diagonal matrix on processor %d is\n", pid);

#if defined(DBL_CMPLX)
    for (i = 0; i < lsize; i++) {
        length = diag_mat->nnzrow[i];
        rja    = diag_mat->pj[i];
        ra     = diag_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "(%d,%d,%d) = (%f, %f)  ", pid, i, rja[j], creal(ra[j]), cimag(ra[j]));
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "Local off-diagonal matrix on processor %d is \n", pid);
    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        rja    = offd_mat->pj[i];
        ra     = offd_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "(%d,%d,%d) = (%f, %f)  ", pid, i, rja[j], creal(ra[j]), cimag(ra[j]));
        }
        fprintf(fp, "\n");
    }
#else
    for (i = 0; i < lsize; i++) {
        length = diag_mat->nnzrow[i];
        rja    = diag_mat->pj[i];
        ra     = diag_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "(%d,%d,%d) = %f  ", pid, i, rja[j], ra[j]);
        }
        fprintf(fp, "\n");
    }

    fprintf(fp, "Local off-diagonal matrix on processor %d is \n", pid);
    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        rja    = offd_mat->pj[i];
        ra     = offd_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "(%d,%d,%d) = %f  ", pid, i, rja[j], ra[j]);
        }
        fprintf(fp, "\n");
    }
#endif
    fprintf(fp, "\n");
    parms_ViewerStoreFP(v, fp);
    parms_CommView(data->mvhandler, v);
    return 0;
}


int parms_MatViewCOO_dvcsr(parms_Mat self, parms_Viewer v)
{
    int i, j, lsize, pid, npro, length;
    int *rja;

    FLOAT *ra;
    parms_dvcsr data;
    parms_vcsr diag_mat, offd_mat;
    parms_Map is;
    FILE *fp;

    is      = self->is;
    lsize   = is->lsize;
    pid     = is->pid;
    npro    = is->npro;
    data    = (parms_dvcsr)self->data;

    parms_ViewerGetFP(v, &fp);

    if (pid == 0) {
        fprintf(fp, "There are %d processors available\n", npro);
        fprintf(fp, "The format of output local equations is:\n");
        fprintf(fp, "local_row_index  local_column_index   value\n");
    }

    fprintf(fp, "\n");

    diag_mat = data->diag_mat;
    offd_mat = data->offd_mat;

#if defined(DBL_CMPLX)
    fprintf(fp, "%d %d %d\n", lsize, lsize, nnzCS(diag_mat) );//getchar();

    for (i = 0; i < lsize; i++) {
        length = diag_mat->nnzrow[i];
        rja    = diag_mat->pj[i];
        ra     = diag_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "%d  %d  %f  %f  \n", i, rja[j], creal(ra[j]), cimag(ra[j]));
        }
    }
    fprintf(fp, "%d %d %d\n", is->ninf, is->ninf, nnzCS(offd_mat) );//getchar();
    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        rja    = offd_mat->pj[i];
        ra     = offd_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "%d  %d  %f  %f  \n", i+lsize-is->ninf, rja[j], creal(ra[j]), cimag(ra[j]));
        }
    }
#else
    fprintf(fp, "%d %d %d\n", lsize, lsize, nnzCS(diag_mat) );//getchar();
    for (i = 0; i < lsize; i++) {
        length = diag_mat->nnzrow[i];
        rja    = diag_mat->pj[i];
        ra     = diag_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "%d  %d  %f  \n", i, rja[j], ra[j]);
        }
    }
    fprintf(fp, "%d %d %d\n", is->ninf, is->ninf, nnzCS(offd_mat) );//getchar();
    for (i = 0; i < is->ninf; i++) {
        length = offd_mat->nnzrow[i];
        rja    = offd_mat->pj[i];
        ra     = offd_mat->pa[i];
        for (j = 0; j < length; j++) {
            fprintf(fp, "%d  %d  %f  \n", i+lsize-is->ninf, rja[j], ra[j]);
        }
    }
#endif
    parms_ViewerStoreFP(v, fp);
    return 0;
}


int parms_MatCreate_dvcsr(parms_Mat self)
{
    parms_dvcsr  data;
    parms_Comm   mvhandler;
    MPI_Comm     comm;

    PARMS_MEMCPY(self->ops, &parms_matops_dvcsr, 1);
    PARMS_NEW(data);
    comm = self->is->comm;
    parms_CommCreate(&mvhandler, comm);
    data->mvhandler = mvhandler;

    self->data = data;
    return 0;
}



//--------------------------------------------------------------------------------------dividing line------------------------------------------------------------------------------------------



/*
typedef struct parms_dvcsr {
  parms_vcsr    diag_mat;
  parms_vcsr    offd_mat;

  parms_bvcsr   b_diag_mat;//vbsptr
  parms_bvcsr   b_offd_mat;//vbsptr
  parms_Comm    mvhandler;

} *parms_dvcsr;*/
/** 
 * Setup the matrix
 *
 * This is the most important function in the package. The local
 * variables are divided into two categories: interior variables (not
 * receiving data from other processors) and interface variables
 * (receiving data from other processors). The interior variables are
 * divided into independent variables (not sending data to other
 * processors) and output variables (sending data to other processors)
 * for *asymmetric* pattern. The independent variables are listed first,
 * followd by output variables, and interface variables.  The member
 * variable schur_start is the number of independent
 * variables. i.e., the start position of "Schur" complement in the
 * local linear system.  communication handler for the matrix-vector
 * product is also created.
 *
 * @param self A matrix object.
 *
 * @return 0 on success.
 */
static int MatSetup_b_dvcsr(parms_Mat self)//when it comes to free of aux_data and diag and offd, i need to go back to this routine and check
{
    parms_dvcsr data;
    parms_bvcsr  b_diag_mat = NULL, b_offd_mat = NULL;
    parms_bvcsr b_aux_data;
    parms_Map  is;
    parms_Comm handler;
    int i, j, *perm, *iperm, mypid, pid, npro, lsize, k1, k2;
    int *rja=NULL, *tja=NULL, nnz, maxnnz, start, end, nrow, nnzb;
    int gnodv, cindex, index, gindex;
    int *ginfvars, *ginfv2p, *ginfv2ps, *odperm, *neword, *clp, *info,
            *disp, *ja;
    int npsendrecv, pos, *procslist=NULL;
    int *ptrvsend, *ptrvrecv=NULL, *cmap=NULL;
    BData *ra=NULL, *tma=NULL, val;//*ra=NULL,FLOAT

    int *lnB = NULL, *newlnB = NULL, *newlnBC = NULL, *ptrvsend_b_nB = NULL;//nB is global dimension of each block, lnB is local; newlnb is new one after local permutation, newlnBC is for the column index of offd

    MPI_Comm comm;

    //output_intvectorpa("is->nBbegin.coo",is->nB,0, is->gsize);getchar();
    /* extract information arrays */
    is            = self->is;
    npro          = is->npro;
    mypid         = is->pid;
    self->issetup = true;
    data          = (parms_dvcsr)self->data;
    b_aux_data    = self->b_aux_data;
    comm          = is->comm;
    handler       = data->mvhandler;
    lsize         = parms_MapGetLocalSize(is);
    perm          = is->perm;
    iperm         = is->iperm;

    //output_intvectorpa("b_aux_data->bsz",b_aux_data->bsz,0, lsize+1);getchar();

    /* free member space */
    if(self->b_aux_data->space)
        PARMS_FREE(self->b_aux_data->space);

    /* allocate memory for the diagonal matrix */
    if(self->isreset)
    {
        if(self->resetpattern == SAME_NONZERO_STRUCTURE)
        {
            /* copy data into parallel matrix struct */
            b_offd_mat = data->b_offd_mat;
            b_diag_mat = data->b_diag_mat;
            /* diagonal part */
            for(i=0; i<b_diag_mat->n; i++)
            {
                PARMS_MEMCPY(b_diag_mat->ba[i],  b_aux_data->ba[i], b_diag_mat->nzcount[i]);
            }
            /* off diagonal part */
            for(i=0; i<b_offd_mat->n; i++)
            {
                index = i+is->nint;
                nrow = b_diag_mat->nzcount[index];
                PARMS_MEMCPY(b_offd_mat->ba[i], &b_aux_data->ba[index][nrow], b_offd_mat->nzcount[i]);
            }
            return 0;
        }
        else
        {
            /* free b_offd_mat for previous data. No need to free b_diag_mat
        * for previous data, since array size (lsize) does not change
        */
            b_diag_mat = data->b_diag_mat;
            b_offd_mat = data->b_offd_mat;
            if(b_offd_mat->n){
                PARMS_FREE(b_offd_mat->nzcount);
                PARMS_FREE(b_offd_mat->ba);
                PARMS_FREE(b_offd_mat->ja);
                //-----------------------------------------------------------------------for block version only
                PARMS_FREE(b_offd_mat->bsz);
                if(b_offd_mat->bszc)
                    PARMS_FREE(b_offd_mat->bszc);
                if(b_offd_mat->D)
                    PARMS_FREE(b_offd_mat->D);
                //-------------------------------------------------------------------------
            }
            PARMS_FREE(b_offd_mat);

            /* free previous communication handler */
            parms_CommFree(&data->mvhandler);
            /* now create new communication handler */
            parms_Comm new_handler;
            parms_CommCreate(&new_handler, comm);
            handler = data->mvhandler = new_handler;
            /* free processor info and table data */
            if(is->vstable) parms_TableFree(&is->vstable);
            if(is->vsend) PARMS_FREE(is->vsend);
        }
    }
    else{
        PARMS_NEW(data->b_diag_mat);
        b_diag_mat = data->b_diag_mat;
        PARMS_NEWARRAY(b_diag_mat->nzcount, lsize);
        PARMS_NEWARRAY(b_diag_mat->ja,     lsize);
        PARMS_NEWARRAY(b_diag_mat->ba,     lsize);
        b_diag_mat->n = lsize;
        //------------------------------------------------------------------------
        b_diag_mat->bsz = NULL;
        b_diag_mat->D = NULL;
        b_diag_mat->bszc = NULL;
        b_diag_mat->lsize = 0;
        //--------------------------------------------------------------------------
    }

    is->nint = lsize-is->ninf;
    is->schur_start = is->nint;
    k1 = 0;
    k2 = is->nint;
    //printf("before perm");
    for (i = 0; i < lsize; i++) {
        if (perm[i] == -1) {
            perm[i] = k1++;
        }
        else {
            perm[i] = k2++;
        }
    }

    /* set up inverse permutation array */
    for (i = 0; i < lsize; i++) {
        k1 = perm[i];
        iperm[k1] = i;
    }
    is->isperm = true;

    /* permute the rows of the local matrix *///permuted b_diag_mat first, then assign to b_aux_data, so permuted both
    for (i = 0; i < lsize; i++) {
        if (b_aux_data->nzcount[i]) {
            PARMS_RESIZE(b_aux_data->ja[i], b_aux_data->nzcount[i]);
            PARMS_RESIZE(b_aux_data->ba[i], b_aux_data->nzcount[i]);
        }
        k1 = perm[i];
        b_diag_mat->ja[k1]     = b_aux_data->ja[i];
        b_diag_mat->ba[k1]     = b_aux_data->ba[i];
        b_diag_mat->nzcount[k1] = b_aux_data->nzcount[i];
    }

    PARMS_MEMCPY(b_aux_data->ja,     b_diag_mat->ja,     lsize);
    PARMS_MEMCPY(b_aux_data->ba,     b_diag_mat->ba,     lsize);
    PARMS_MEMCPY(b_aux_data->nzcount, b_diag_mat->nzcount, lsize);

    /* permute the column indices of the matrix which are corresponding
     to the interior variables *///only b_aux_data
    for (i = 0; i < is->nint; i++) {
        rja = b_aux_data->ja[i];
        nnz = b_aux_data->nzcount[i];
        ra  = b_aux_data->ba[i];
        for (j = 0; j < nnz; j++) {
            rja[j] = perm[rja[j]];	/* in C style */
        }
    }

    /* permute the local column indices of the matrix which are
     corresponding to the interface variables *///also b_aux_data
    maxnnz = 0;
    for (i = is->nint; i < lsize; i++) {
        rja = b_aux_data->ja[i];
        nnz = b_aux_data->nzcount[i];
        if (maxnnz < nnz) {
            maxnnz =  nnz;
        }
        b_diag_mat->nzcount[i] = 0;
        for (j = 0; j < nnz; j++) {
            if (rja[j] >= 0) {
                b_diag_mat->nzcount[i]++;
                rja[j] = perm[rja[j]]; /* in C style */
            }
        }
    }

    if (maxnnz) {
        PARMS_NEWARRAY(tja, maxnnz);
        PARMS_NEWARRAY(tma, maxnnz);
    }

    /* allocate memory for the off-diagonal matrix *///??
    PARMS_NEW(data->b_offd_mat);
    b_offd_mat = data->b_offd_mat;
    b_offd_mat->n = is->ninf;
    if (is->ninf) {
        PARMS_NEWARRAY(b_offd_mat->nzcount, is->ninf);
        PARMS_NEWARRAY(b_offd_mat->ja,     is->ninf);
        PARMS_NEWARRAY(b_offd_mat->ba,     is->ninf);
        //---------------------------------------------------------------------
        b_offd_mat->bsz = NULL;
        b_offd_mat->D = NULL;
        b_offd_mat->bszc = NULL;
        b_offd_mat->lsize = 0;

        //---------------------------------------------------------------------
    }
    handler->nodv = parms_TableGetSize(is->table) - lsize;
    if (handler->nodv) {
        PARMS_NEWARRAY(handler->odvlist, handler->nodv);
    }

    /* list column indices of diagonal part first */
    for (i = is->nint; i < lsize; i++) {
        rja    = b_aux_data->ja[i];
        ra     = b_aux_data->ba[i];
        nnz    = b_aux_data->nzcount[i];
        index  = i-is->nint;
        k1     = 0;
        k2     = b_diag_mat->nzcount[i];
        nrow   = b_diag_mat->nzcount[i];
        for (j = 0; j < nnz; j++) {
            if (rja[j] >= 0) {
                tja[k1]   = rja[j];
                tma[k1++] = ra[j];
            }
            else {
                tja[k2]   = -rja[j]-1;
                clp = parms_MapGlobalToLocal(is, tja[k2]);
                if (clp != NULL) {
                    handler->odvlist[*clp-lsize] = tja[k2];
                }
                tma[k2++] = ra[j];
            }
        }
        PARMS_MEMCPY(rja, tja, nnz);//(dse, sou, nnz)
        PARMS_MEMCPY(ra,  tma, nnz);
        b_offd_mat->ja[index]     = &rja[nrow];
        b_offd_mat->ba[index]     = &ra[nrow];
        b_offd_mat->nzcount[index] = nnz-nrow;
    }
    if (maxnnz) {
        PARMS_FREE(tja);
        PARMS_FREE(tma);
    }

    /* sort the column indices of the diagonal matrix in ascending
     order */
    for (i = 0; i < lsize; i++) {
        nnz = b_diag_mat->nzcount[i];
        rja = b_diag_mat->ja[i];
        ra  = b_diag_mat->ba[i];
        for (j = 0; j < nnz-1; j++) {
            for (k1 = j+1; k1 < nnz; k1++) {
                if (rja[j] > rja[k1]) {
                    index = rja[j];
                    rja[j] = rja[k1];
                    rja[k1] = index;

                    val = ra[j];
                    ra[j] = ra[k1];
                    ra[k1] = val;
                }
            }
        }
    }

    PARMS_NEWARRAY(info, npro);
    PARMS_NEWARRAY(disp, npro+1);

    MPI_Allgather(&handler->nodv, 1, MPI_INT, info, 1, MPI_INT, comm);
    disp[0] = 0;
    for (i = 0; i < npro; i++) {
        disp[i+1] = disp[i] + info[i];
    }

    /* total number of interface variables */
    gnodv = disp[npro];
    if (gnodv) {
        PARMS_NEWARRAY(ginfvars,  gnodv);
        PARMS_NEWARRAY0(ginfv2p,  gnodv);
        PARMS_NEWARRAY(ginfv2ps,  gnodv);
        MPI_Allgatherv(handler->odvlist, handler->nodv, MPI_INT, ginfvars,
                       info, disp, MPI_INT, comm);

        /* variables sent to other processors */
        for (i = 0; i < npro; i++) {
            if (i != mypid) {
                for (j = disp[i]; j < disp[i+1]; j++) {
                    clp = parms_MapGlobalToLocal(is, ginfvars[j]);
                    if (clp != NULL && *clp < lsize) {
                        ginfv2p[j] = mypid;
                    }
                }
            }
        }

        MPI_Allreduce(ginfv2p, ginfv2ps, gnodv, MPI_INT, MPI_SUM, comm);
        PARMS_FREE(ginfv2p);
        PARMS_NEWARRAY(handler->procs_recv, npro);
        for (i = 0; i < npro; i++) {
            info[i] = 0;
        }
        /* processor IDs from which processor mypid receives data */
        handler->nprecv = 0;
        for (i = disp[mypid]; i < disp[mypid+1]; i++) {
            pid = ginfv2ps[i];
            if (info[pid] == 0) {
                handler->procs_recv[handler->nprecv++] = pid;
            }
            info[pid]++;
        }

        if (handler->nprecv) {
            PARMS_RESIZE(handler->procs_recv, handler->nprecv);

            /* sort neighbouring processors received data from in ascending order */
            for (i = 0; i < handler->nprecv-1; i++) {
                for (j = i+1; j < handler->nprecv; j++) {
                    if (handler->procs_recv[i] > handler->procs_recv[j]) {
                        pid = handler->procs_recv[i];
                        handler->procs_recv[i] = handler->procs_recv[j];
                        handler->procs_recv[j] = pid;
                    }
                }
            }
        }
        else {
            PARMS_FREE(handler->procs_recv);
        }

        PARMS_NEWARRAY(handler->ptrvrecv, handler->nprecv+1);
        handler->ptrvrecv[0] = 0;
        for (i = 0; i < handler->nprecv; i++) {
            pid = handler->procs_recv[i];
            handler->ptrvrecv[i+1] = handler->ptrvrecv[i] + info[pid];
            info[pid] = i;
        }

        if (handler->nodv) {
            PARMS_NEWARRAY(odperm, handler->nodv);
            //PARMS_NEWARRAY(handler->buf_recv, handler->nodv);

            /* reorder external variables: the data received from the same
     processor are listed consecutively, processor by processor */
            j = disp[mypid];
            for (i = 0; i < handler->nodv; i++) {
                pid = ginfv2ps[i+j];
                index = info[pid];
                odperm[i] = handler->ptrvrecv[index]++;
            }
            for (i = handler->nprecv-1; i > 0; i--) {
                handler->ptrvrecv[i] = handler->ptrvrecv[i-1];
            }
            handler->ptrvrecv[0] = 0;

            /* reorder the external interface variables */
            PARMS_NEWARRAY(neword, handler->nodv);
            for (i = 0; i < handler->nodv; i++) {
                index = odperm[i];
                gindex = handler->odvlist[i];
                parms_TablePut(is->table, gindex, index+lsize);// insert the pair (key, value) into the table.
                neword[index] = gindex;
            }

            PARMS_FREE(handler->odvlist);

            PARMS_FREE(odperm);
            handler->odvlist = neword;
        }
        /* change the global column indices of entries in external
       matrix to the local ones */
        //outputcsmatpa(offd_mat,"offd_datainsetupbefore_local",1);

        for (i = 0; i < b_offd_mat->n; i++) {
            nnz = b_offd_mat->nzcount[i];
            ja  = b_offd_mat->ja[i];
            for (j = 0; j < nnz; j++) {
                cindex = ja[j];
                clp  = parms_MapGlobalToLocal(is, cindex);
                ja[j]  = *clp;

                //if (pid == 0 ) printf("row = %d, cindex = %d, clp = %d \n",i, cindex,*clp);//getchar();

            }
        }
        PARMS_FREE(ginfv2ps);
        PARMS_FREE(ginfvars);
    }
    else {
        handler->nprecv = 0;
        PARMS_NEWARRAY(handler->ptrvrecv, 1);
        handler->ptrvrecv[0] = 0;
    }

    /* set up send information */
    MPI_Allgather(&handler->nprecv, 1, MPI_INT, info, 1, MPI_INT, comm);
    disp[0] = 0;
    for (i = 0; i < npro; i++) {
        disp[i+1] = disp[i] + info[i];
    }

    if (disp[npro]) {
        PARMS_NEWARRAY(procslist, disp[npro]);
        MPI_Allgatherv(handler->procs_recv, handler->nprecv, MPI_INT,
                       procslist, info, disp, MPI_INT, comm);
    }

    handler->npsend = 0;
    PARMS_NEWARRAY(handler->procs_send, npro);
    for (i = 0; i < npro; i++) {
        if (i != mypid) {
            for (j = disp[i]; j < disp[i+1]; j++) {
                if (procslist[j] == mypid) {
                    handler->procs_send[handler->npsend++] = i;
                    break;
                }
            }
        }
    }

    if (disp[npro]) {
        PARMS_FREE(procslist);
    }

    if (handler->npsend) {
        PARMS_RESIZE(handler->procs_send, handler->npsend);
    }
    else {
        PARMS_FREE(handler->procs_send);
    }

    PARMS_NEWARRAY(handler->ptrvsend, handler->npsend+1);
    handler->ptrvsend[0] = 0;
    ptrvsend = handler->ptrvsend;
    ptrvrecv = handler->ptrvrecv;
    npsendrecv = handler->npsend + handler->nprecv;
    if (npsendrecv) {
        PARMS_NEWARRAY(handler->status_send, npsendrecv);
        PARMS_NEWARRAY(handler->req_send,    npsendrecv);
        handler->status_recv = handler->status_send + handler->npsend;
        handler->req_recv    = handler->req_send + handler->npsend;
        for (i = 0; i < handler->npsend; i++) {
            MPI_Irecv(&disp[i], 1, MPI_INT, handler->procs_send[i],
                      100, comm, &handler->req_send[i]);
        }

        for (i = 0; i < handler->nprecv; i++) {
            info[i] = ptrvrecv[i+1] - ptrvrecv[i];
            MPI_Isend(&info[i], 1, MPI_INT, handler->procs_recv[i],
                      100, comm, &handler->req_recv[i]);
        }

        MPI_Waitall(handler->nprecv, handler->req_recv,
                    handler->status_recv);
        MPI_Waitall(handler->npsend, handler->req_send,
                    handler->status_send);

        for (i = 0; i < handler->npsend; i++) {
            ptrvsend[i+1] = ptrvsend[i] + disp[i];
        }

        nnz = ptrvsend[handler->npsend];
        if (nnz) {
            //PARMS_NEWARRAY(handler->buf_send,   nnz);
            PARMS_NEWARRAY(handler->vlist_send, nnz);
        }

        handler->mdata_send = 0;
        if (handler->npsend) {
            for (i = 0; i < handler->npsend; i++) {
                nnz = ptrvsend[i+1] - ptrvsend[i];
                if (handler->mdata_send < nnz) {
                    handler->mdata_send = nnz;
                }
            }
        }

        /* handler->vlist_send stores the variables to be sent to adjacent
       processors in global labeling */
        for (i = 0; i < handler->npsend; i++) {
            pos = ptrvsend[i];
            nnz = ptrvsend[i+1] - ptrvsend[i];
            MPI_Irecv(&handler->vlist_send[pos], nnz, MPI_INT,
                      handler->procs_send[i], 200, comm,
                      &handler->req_send[i]);
        }

        for (i = 0; i < handler->nprecv; i++) {
            pos = ptrvrecv[i];
            nnz = ptrvrecv[i+1] - ptrvrecv[i];
            MPI_Isend(&handler->odvlist[pos], nnz, MPI_INT,
                      handler->procs_recv[i], 200, comm,
                      &handler->req_recv[i]);
        }

        MPI_Waitall(handler->nprecv, handler->req_recv,
                    handler->status_recv);
        MPI_Waitall(handler->npsend, handler->req_send,
                    handler->status_send);

        /* change global labelings  to local ones */
        is->ninf_send = 0;
        parms_TableCreate(&is->vstable, NULL, is->ninf);
        if (is->isperm) {
            for (i = 0; i < handler->npsend; i++) {
                for (j = ptrvsend[i]; j < ptrvsend[i+1]; j++) {
                    cindex = handler->vlist_send[j];
                    clp  = parms_MapGlobalToLocal(is, cindex);
                    index  = perm[*clp];
                    handler->vlist_send[j] = index;
                    clp = parms_TableGet(is->vstable, index);
                    if (clp == NULL) {
                        parms_TablePut(is->vstable, index, is->ninf_send++);
                    }
                }
            }
        }

        /* vsend is an array to store variables sent to other processors */
        if (is->ninf_send) {
            PARMS_NEWARRAY(is->vsend, is->ninf_send);
            if (is->nint) {
                PARMS_NEWARRAY(cmap, is->nint);
            }
            for (i = 0; i < is->nint; i++) {
                cmap[i] = -1;
            }

            for (i = 0; i < handler->npsend; i++) {
                for (j = ptrvsend[i]; j < ptrvsend[i+1]; j++) {
                    index = handler->vlist_send[j];
                    clp = parms_TableGet(is->vstable, index);
                    is->vsend[*clp] = index;

                }
            }
            for (i = 0; i < is->ninf_send; i++) {
                if (is->vsend[i] < is->nint) {
                    cmap[is->vsend[i]] = --is->schur_start;
                }
            }

            /* reorder the local variables in the following order:
       * 1. independent variables which are uncoupled with variables in
       *    adjacent processors.
       * 2. variables sent to adjacent processors but the
       *    corresponding rows are decoupled with variables in
       *    the adjacent processors. (asymmetric pattern)
       * 3. interface variables - coupled with variables in adjacent
       *    processors in the sense of receiving data from
       *    the adjacent processors.
       */
            if (is->schur_start != is->nint) {
                k1 = 0;
                for (i = 0; i < is->nint; i++) {
                    if (cmap[i] == -1) {
                        cmap[i] = k1++;
                    }
                }
                /* permute column indices of the local matrix */
                for (i = 0; i < lsize; i++) {
                    nnz = b_diag_mat->nzcount[i];
                    ja = b_diag_mat->ja[i];
                    for (j = 0; j < nnz; j++) {
                        if (ja[j] < is->nint) {
                            ja[j] = cmap[ja[j]];
                        }
                    }
                }

                /* permute rows which are corresponding to interior variables */
                for (i = 0; i < is->nint; i++) {
                    k1 = cmap[i];
                    b_diag_mat->nzcount[k1] = b_aux_data->nzcount[i];
                    b_diag_mat->ja[k1]     = b_aux_data->ja[i];
                    b_diag_mat->ba[k1]     = b_aux_data->ba[i];
                }
                PARMS_MEMCPY(b_aux_data->nzcount, b_diag_mat->nzcount, is->nint);
                PARMS_MEMCPY(b_aux_data->ja,     b_diag_mat->ja,     is->nint);
                PARMS_MEMCPY(b_aux_data->ba,     b_diag_mat->ba,     is->nint);

                /* update perm and iperm */
                for (i = 0; i < lsize; i++) {
                    if (perm[i] < is->nint) {
                        perm[i] = cmap[perm[i]];
                    }
                }

                for (i = 0; i < lsize; i++) {
                    k1 = perm[i];
                    iperm[k1] = i;
                }

                /* update vlist_send */
                for (i = 0; i < handler->npsend; i++) {
                    start = handler->ptrvsend[i];
                    end   = handler->ptrvsend[i+1];
                    for (j = start; j < end; j++) {
                        if (handler->vlist_send[j] < is->nint) {
                            handler->vlist_send[j] = cmap[handler->vlist_send[j]];
                        }
                    }
                }

                /* update vsend */
                for (i = 0; i < is->ninf_send; i++) {
                    if (is->vsend[i] < is->nint) {
                        is->vsend[i] = cmap[is->vsend[i]];
                        parms_TablePut(is->vstable, is->vsend[i], i);
                    }
                }
            }
            if (is->nint) {
                PARMS_FREE(cmap);
            }

        }

    }

    //printf("before bsz permutation, lsize = %d, b_aux_data->bsz[0] = %d, b_aux_data->bsz[1376] = %d\n",lsize, b_aux_data->bsz[0],b_aux_data->bsz[1376]);getchar();
    //for (i=0; i<lsize+1; i++) {
    //printf("b_aux_data->bsz[%d] = %d, mypid = %d\n",i, b_aux_data->bsz[i],mypid);//getchar();
    //}
    //output_intvectorpa("b_aux_data->bsz",b_aux_data->bsz,0, lsize+1);//getchar();

    /*this beneath is to contruct the block structure of b_diag, b_offd, b_aux_data*/
    /*and also set up the data communication handler members*/

    /*first generate the b_diag_mat->bsz by permuting the old one*/

    PARMS_NEWARRAY(b_diag_mat->bsz,     lsize+1);
    b_diag_mat->bsz[0] = 0;

    PARMS_NEWARRAY(lnB,lsize);//allocate memory to local block structure array, only for parallel case, vbsptr only
    PARMS_NEWARRAY(newlnB,lsize);

    for (i=0; i<lsize; i++) {
        lnB[i] = B_DIM(b_aux_data->bsz,i);
        newlnB[perm[i]] = lnB[i];
    }

    //printf()
    //output_intvectorpa("lperm",perm,0, lsize);//getchar();
    //output_intvectorpa("lnB",lnB,0, lsize);//getchar();
    //output_intvectorpa("newlnB",newlnB,0, lsize);//getchar();
    for (i = 0; i < lsize; i++) {
        b_diag_mat->bsz[i+1] = b_diag_mat->bsz[i] + newlnB[i];
    }
    is->lbsz_after = b_diag_mat->bsz;

    /*update the b_aux_data->bsz*/
    PARMS_MEMCPY( b_aux_data->bsz,b_diag_mat->bsz, lsize+1);

    /* output_intvectorpa("b_diag_mat->bsz",b_diag_mat->bsz,0, b_diag_mat->n+1);//getchar(); */
    /* output_intvectorpa("b_aux_data->bsz",b_aux_data->bsz,0, b_aux_data->n+1);//getchar(); */
    /* MPI_Barrier(MPI_COMM_WORLD); */
    /* exit(1); */
    //-------------------------------------------------------------------------------------------
    PARMS_NEWARRAY(is->pperm,is->llsize);
    permuten(perm, is->pperm, b_diag_mat->bsz, lsize);  //permuten(map->perm, pperm, map->lbsz, map->lsize);//int permuten(int *iwork, int *iworkn, int *bsz, int n)
    //  for (i = 0; i < lsize+1; i++) {
    //    is->lbsz[i] = b_diag_mat->bsz[i];
    //  }//store local block structure in map struct
    //printf(" handler->nodv = %d, boffd = %d, lsize = %d\n",handler->nodv, b_offd_mat, lsize);getchar();
    //printf("  vlist_send= %d, mypid = %d\n",handler->vlist_send, mypid);

    /*set up some handler member */
    /*first is the members for sending data out*/

    if(handler->vlist_send){

        PARMS_NEWARRAY(handler->ptrvsend_b, handler->npsend+1);
        handler->ptrvsend_b[0] = 0;
        //ptrvsend_b = handler->ptrvsend_b;
        PARMS_NEWARRAY0(ptrvsend_b_nB, handler->npsend);

        for (i = 0; i < handler->npsend; i++) {
            //fprintf(fp, "The data in local indices sent to processor %d\n", handler->procs_send[i]);
            for (j = handler->ptrvsend[i]; j < handler->ptrvsend[i+1]; j++) {
                ptrvsend_b_nB[i] += newlnB[handler->vlist_send[j]];
                //fprintf(fp, "vlist_send[%d] = %d\n", j, handler->vlist_send[j]);

            }
            //printf(" ptrvsend_b_nB[%d] = %d, handler->npsend = %d, mypid = %d\n",i,ptrvsend_b_nB[i], handler->npsend, mypid);
            handler->ptrvsend_b[i+1] = handler->ptrvsend_b[i] + ptrvsend_b_nB[i];
        }



        nnzb = handler->ptrvsend_b[handler->npsend];//point version total number of variabls need to be sent to other processors
        //nnz = handler->ptrvsend_b[handler->npsend];
        //if (nnz) {
        PARMS_NEWARRAY(handler->buf_send_b,   nnzb);// to store the value of the data
        //PARMS_NEWARRAY(handler->vlist_send, nnzb);
        //}
        PARMS_NEWARRAY(handler->vlist_send_b,  nnzb);// to store the index of the data
        nnz = handler->ptrvsend[handler->npsend];//block version total number of variabls need to be sent to other processors

        permuten(handler->vlist_send, handler->vlist_send_b, b_diag_mat->bsz, nnz);//int permuten(int *iwork, int *iworkn, int *bsz, int n)
        //for (i = 0; i < nnz; i++) {
        //handler->lbsz[i] = b_diag_mat->bsz[i];
        //}//store local block structure in map struct

        //for (i = 0; i < lsize+1; i++) {
        //handler->lbsz[i] = b_diag_mat->bsz[i];
        //}//store local block structure in map struct
        //output_intvectorpa("handler->lbsz",handler->lbsz,0, lsize+1);//getchar();
        //vlist_send_b

        //output_intvectorpa("handler->ptrvsend_b",handler->ptrvsend_b,0, handler->npsend+1);//getchar();
        //output_intvectorpa("handler->ptrvsend",handler->ptrvsend,0, handler->npsend+1);//getchar();
        //output_intvectorpa("handler->vlist_send",handler->vlist_send,0, nnz);//getchar();
        //output_intvectorpa("handler->vlist_send_b",handler->vlist_send_b,0, nnzb);//getchar();
        PARMS_FREE(ptrvsend_b_nB);

    }
    //MPI_Barrier(MPI_COMM_WORLD);
    //------------------------------------------------------------------------------


    //for (i = 0; i < ptrvsend[handler->npsend+1]; i++) {
    //sum += newlnB[vlist_send[i]];//handler->lbsz[i] = b_diag_mat->bsz[i];
    //}//store local block structure in map struct
    //printf(" sum = %d\n",sum);



    //handler->ptrvrecv_b[0] = 0;
    //ptrvrecv_b = handler->ptrvrecv_b;

    /*set up some communication handler member */
    /*and also insert value into b_offd_mat->bsz, b_offd_mat->bszc and b_aux_data->bszc*/
    /*actually, b_aux_data = b_offd_mat append to b_diag_mat*/
    /*set up the members for receiving data from other processors*/

    if (handler->nodv) {
        b_offd_mat->nc = handler->nodv;
        b_offd_mat->lsize = lsize;//this is block column index offset

        //if(handler->nodv) {
        //printf(" b_offd_mat->nc = %d\n",b_offd_mat->nc);
        //}


        PARMS_NEWARRAY(newlnBC,handler->nodv);//specially for offd, bszc
        PARMS_NEWARRAY(b_offd_mat->bsz,     is->ninf+1);
        b_offd_mat->bsz[0] = 0;
        PARMS_NEWARRAY(b_offd_mat->bszc,     handler->nodv+1);
        b_offd_mat->bszc[0] = 0;

        for (i = 0; i < is->ninf; i++) {
            b_offd_mat->bsz[i+1] = b_offd_mat->bsz[i] + newlnB[is->nint+i];
        }

        for (i = 0; i < handler->nodv; i++) {
            newlnBC[i] = is->nB[handler->odvlist[i]];//gindex = handler->odvlist[cindex-lsize];
        }

        //output_intvectorpa("newlnBC",newlnBC,0, b_offd_mat->nc);//getchar();

        int *bszc;
        b_aux_data->nc = lsize + handler->nodv;

        //printf(" handler->nodv = %d\n",handler->nodv);

        PARMS_NEWARRAY(b_aux_data->bszc, lsize+handler->nodv+1); //the length is lsize+handler->nodv+1
        PARMS_MEMCPY(b_aux_data->bszc, b_aux_data->bsz, lsize+1);//update the b_aux_data->bsz
        bszc = b_aux_data->bszc;

        for (i = 0; i < handler->nodv; i++) {
            b_offd_mat->bszc[i+1] = b_offd_mat->bszc[i] + newlnBC[i];//fill in b_offd_mat->bszc
            bszc[i+lsize+1] = b_offd_mat->bszc[i+1] + bszc[lsize];//fill in the second half of b_aux_data->bszc based on b_offd_mat->bszc
        }


        /* for (i = 0; i < handler->nodv; i++) { */
        /*   b_offd_mat->bszc[i+1] = b_offd_mat->bszc[i] + newlnBC[i]; */
        /* } */
        //---------------------------------------------------------------------
        PARMS_NEWARRAY(handler->ptrvrecv_b, handler->nprecv+1);
        //handler->ptrvrecv_b[0] = 0;
        for (i = 0; i < handler->nprecv+1; i++) {
            handler->ptrvrecv_b[i] = b_offd_mat->bszc[ptrvrecv[i]];
            //printf(" ptrvrecv_b[%d] = %d, handler->nodv = %d, mypid = %d\n",i,handler->ptrvrecv_b[i], handler->nodv, mypid);
        }
        //MPI_Barrier(MPI_COMM_WORLD);

        nnz = handler->ptrvrecv_b[handler->nprecv];
        //if (nnz) {
        PARMS_NEWARRAY(handler->buf_recv_b,   nnz);
        //PARMS_NEWARRAY(handler->vlist_send, nnz);
        //}

        //output_intvectorpa("handler->ptrvrecv_b",handler->ptrvrecv_b,0, handler->nprecv+1);//getchar();
        //output_intvectorpa("handler->ptrvrecv",handler->ptrvrecv,0, handler->nprecv+1);//getchar();

        //-----------------------------------------------------------------------
        //PARMS_NEWARRAY0(ptrvrecv_b_nB, handler->nprecv);
        //output_intvectorpa("b_offd_mat->bsz",b_offd_mat->bsz,0, b_offd_mat->n+1);//getchar();
        //output_intvectorpa("b_offd_mat->bszc",b_offd_mat->bszc,0, b_offd_mat->nc+1);//getchar();
        /* outputvbmatpa(b_offd_mat,"b_offd_data",1); */
        //output_intvectorpa("columnglobal",handler->odvlist,0, handler->nodv);
        //MPI_Barrier(MPI_COMM_WORLD);
        //exit(1);
    }
    /* else { */
    /*   b_aux_data->nc = b_aux_data->n; */
    /*   b_aux_data->bszc = b_aux_data->bsz; */
    /* } */
    //output_intvectorpa("localpermafter",perm,0, lsize);//getchar();
    //output_intvectorpa("b_diag_mat->bsz",b_diag_mat->bsz,0, b_diag_mat->n+1);//getchar();
    //output_intvectorpa("is->pperm",is->pperm,0, is->llsize);//getchar();
    free(lnB);lnB = NULL;
    free(newlnB);newlnB = NULL;
    if (newlnBC) {
        free(newlnBC);newlnBC = NULL;
    }

    //PARMS_FREE(ptrvrecv_b_nB);
    //for (i = 0; i < handler->nodv; i++) {
    //b_offd_mat->bszc[i+1] = b_offd_mat->bszc[i] + newlnB[is->nint+i];
    //}
    //outputvbmatpa(b_diag_mat,"b_diag_data",1);
    //outputvbmatpa(b_offd_mat,"b_offd_data",1);
    //outputvbmatpa(b_aux_data,"b_aux_data",1);
    //parms_VecPerm(x, is);


    //parms_VecPerm(x, is);
    //parms_VecInvPerm(x, is);
    PARMS_FREE(info);
    PARMS_FREE(disp);

    return 0;
}

static int MatVec_b_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y) 
{
    int lsize, i, j, *ja, length, dim, sz, nBs, nBsj, col, inc = 1, row;//, llsize;
    int *bsz, *bszc;
    parms_dvcsr data;
    parms_Comm  handler;
    parms_bvcsr b_diag_mat, b_offd_mat;
    parms_Map is;
    FLOAT *offsetptr;//, s;//*ba,
    BData *ba;
    FLOAT one=1.0;

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
//    llsize   = is->llsize;

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
            GGEMV ("n", dim, sz,one, ba[j],dim,&x[nBsj],inc,one,&y[nBs],inc);//s    += ba[j] * x[index];
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
            GGEMV ("n", dim, sz,one, ba[j],dim,&offsetptr[nBsj],inc,one,&y[nBs],inc);//s    += ba[j] * x[index];

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

static int MatMVPY_b_dvcsr(parms_Mat self, FLOAT alpha, FLOAT *x,
                           FLOAT beta, FLOAT *y, FLOAT *z)// not done yet
{

    return 0;
}

static int MatVecSchur_b_dvcsr(parms_Mat self, FLOAT *x, FLOAT *y, int
                               pos)
{
    parms_dvcsr data;
    parms_Comm  handler;
    parms_Map   is;
    //parms_vcsr  offd_mat;
    parms_bvcsr b_diag_mat, b_offd_mat;
    FLOAT       *offsetptr;//, s;//*pa,
    int i, j, length, *ja, lsize, *bsz, *bszc;
    int dim, sz, nBs, nBsj, col, inc = 1, row;//, llsize;
    BData *ba;
    FLOAT one=1.0;

    is       = self->is;
    lsize    = is->lsize;
//    llsize   = is->llsize;
    data     = (parms_dvcsr)self->data;
    handler  = data->mvhandler;
    b_diag_mat = data->b_diag_mat;
    b_offd_mat = data->b_offd_mat;

    if ((x == NULL) || (y == NULL)) {
        return 0;
    }


    bsz = b_diag_mat->bsz;
    bszc = b_offd_mat->bszc;

    //printf("pos = %d\n", pos);
    //pos = bsz[pos];
    //printf("pos = %d\n", pos);

    parms_CommDataBegin_b(handler, x, pos);
    for (i = pos; i < bsz[is->nint]; i++) {//for (i = pos; i < is->nint; i++) {
        y[i-pos] = 0.0;
    }
    parms_CommDataEnd(handler);
    offsetptr = handler->buf_recv_b;
    for (i = 0; i < is->ninf; i++) {

        row = i+is->nint;
        nBs = bsz[row];
        dim = B_DIM(bsz,row);

        for( j = 0; j < dim; j++ )
            y[nBs+j-pos] = 0.0;//s      = 0.0; do something here

        length = b_offd_mat->nzcount[i];
        ja     = b_offd_mat->ja[i];
        ba     = b_offd_mat->ba[i];

        for (j = 0; j < length; j++) {
            col = ja[j]-lsize;
            nBsj = bszc[col];
            //index = ja[j];
            sz = B_DIM(bszc,col);
            //-------------------- operation:  y = Block*x + y
            GGEMV ("n", dim, sz,one, ba[j],dim,&offsetptr[nBsj],inc,one,&y[nBs-pos],inc);//s    += ba[j] * x[index];      //s   += ba[j] * offsetptr[ja[j]];
        }
        //y[i+is->nint-pos] = s;
    }
    return 0;
}

static int MatGetDiag_b_dvcsr(parms_Mat self, void **mat)
{
    parms_dvcsr data;
    parms_bvcsr  b_diag_mat;

    // printf("in MatGetDiag_b_dvcsr");

    data = (parms_dvcsr)self->data;
    b_diag_mat = data->b_diag_mat;
    //if (self->ilutype == PCVBARMS) {//??
    if ((self->ilutype == PCVBARMS || self->ilutype == PCVBARMSOLD) && self->isin_b_schur == false) {
        parms_bvcsr diag;
        int i, j, nnz;
        int dimR, dimC, col, bnz;


        //printf("in MatGetDiag_dvcsr a new diag");

        PARMS_NEW(diag);
        diag->n = b_diag_mat->n;
        diag->bszc = NULL;
        diag->lsize = 0;
        diag->D = NULL;
        PARMS_NEWARRAY(diag->bsz, diag->n+1);
        diag->bsz[0] = 0;
        PARMS_NEWARRAY(diag->nzcount, diag->n);
        PARMS_NEWARRAY(diag->ja,     diag->n);
        PARMS_NEWARRAY(diag->ba,     diag->n);
        for (i = 0; i < diag->n; i++) {
            nnz = diag->nzcount[i] = b_diag_mat->nzcount[i];
            diag->bsz[i+1] = b_diag_mat->bsz[i+1];

            if (nnz) {
                //printf("getdiag");
                dimR = B_DIM(b_diag_mat->bsz,i);
                PARMS_NEWARRAY(diag->ja[i], nnz);
                PARMS_NEWARRAY(diag->ba[i], nnz);
                PARMS_MEMCPY(diag->ja[i], b_diag_mat->ja[i], nnz);
                for (j = 0; j < nnz; j++){
                    col = diag->ja[i][j];
                    dimC = B_DIM(b_diag_mat->bsz,col);
                    bnz = dimC*dimR;
                    if (bnz) {
                        PARMS_NEWARRAY(diag->ba[i][j], bnz);
                        PARMS_MEMCPY(diag->ba[i][j], b_diag_mat->ba[i][j], bnz);
                    }
                }
                //PARMS_MEMCPY(diag->ba[i], b_diag_mat->ba[i], nnz);//??
            }
        }
        *mat = diag;
        return 0;
    }
    //printf("in MatGetDiag_dvcsr use old  diag");

    *mat = b_diag_mat;
    return 0;
}


static int MatGetSubMat_b_dvcsr(parms_Mat self, void **mat)
{
    *mat = self->b_aux_data;
    return 0;
}

static int MatExtend_b_dvcsr(parms_Mat self, parms_Comm handler, int
                             start, void *mat, int *nn, void
                             **ext_mat)
{
    parms_bvcsr lmat, new_mat;
    parms_Map  is;
    MPI_Request *sreq=NULL, *rreq=NULL;
    MPI_Status *sstatus=NULL, *rstatus=NULL;
    MPI_Comm comm;
    // FLOAT **a_send=NULL, *ba_rbuf=NULL, *ba_sbuf=NULL, *rowa=NULL;
    BData **a_send=NULL, *ba_rbuf=NULL, *ba_sbuf=NULL, *rowa=NULL;
    int **ja_send=NULL;
    int *ia_send=NULL, *pia_send=NULL, *pia_recv=NULL, *ja_rbuf=NULL, *ja_sbuf=NULL;
    int *pindex, *rowj, *rbuf=NULL, *sbuf=NULL;
    int lsize, jcol, i, j, index, nnz, nodv, numsend, numrecv;
    int n, pos, lindex, tag, length, k, m;
    int dimR, dimC, *bsz, *bszc;

    int *pia_send_b, *pia_recv_b, jj, indexnB, knz, blocksz, nBsum, kk;
    FLOAT *ba_rbuf_b = NULL, *ba_sbuf_b = NULL;// The array to store the unpacked data from ba_rbuf and ba_sbuf
    int *ba_rbuf_nB = NULL, *ba_sbuf_nB = NULL;// The array to store the length of each block in receiving and sending buffer
    BData blockbuffer;

    lmat = (parms_bvcsr)mat;
    PARMS_NEW(new_mat);
    is = self->is;
    lsize = is->lsize;

    MPI_Comm_dup(handler->comm, &comm);

    bsz = lmat->bsz;
    if (lmat->bszc)
        bszc = lmat->bszc;
    else bszc = bsz;

    //printf("is->ninf_send = %d\n", is->ninf_send);

    int myid;
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    /* printf("is->nint = %d, pid = %d \n",is->nint, myid); */
    /* printf("is->ninf = %d, pid = %d \n",is->ninf, myid); */
    /* printf("is->ninf_send = %d, pid = %d \n",is->ninf_send, myid); */
    //printf("The number of stored entries is %d in local block form matrix, pid = %d \n",memVBMat(Amat), myid);


    if (is->ninf_send) {
        PARMS_NEWARRAY(a_send,  is->ninf_send);//2D pointer to values
        PARMS_NEWARRAY(ja_send, is->ninf_send);//2D pointer to index
        PARMS_NEWARRAY(ia_send, is->ninf_send);//array to store nnz of each row

        /* copy the matrix corresponding to the interface variables to
       working arrays. the column indices are stored in global numbering */
        for (i = 0; i < is->ninf_send; i++) {
            index = is->vsend[i] - start;//the index of variables to send out
            //dimR = B_DIM(bsz, index);
            nnz = ia_send[i] = lmat->nzcount[index];
            PARMS_NEWARRAY(a_send[i],  nnz);
            PARMS_NEWARRAY(ja_send[i], nnz);
            nnz = lmat->nzcount[index];
            for (j = 0; j < nnz; j++) {
                a_send[i][j] = lmat->ba[index][j];
                jcol = lmat->ja[index][j];
                /* dimC = B_DIM(bszc, jcol); //????? */
                /* //~ printf("dimc %d jcol %d n %d i %d\n", dimC, jcol, bsz[lmat->n], index); */
                /* PARMS_NEWARRAY(a_send[i][j], dimC*dimR); */
                /* copyBData(dimC, dimR, a_send[i][j], lmat->ba[index][j], 0); */

                if (jcol < lsize) {
                    if (is->isperm) {
                        jcol = is->iperm[jcol];
                        ja_send[i][j] = is->lvars[jcol];
                    }
                    else {
                        ja_send[i][j] = is->lvars[jcol];
                    }
                }
                else {
                    jcol = handler->odvlist[jcol-lsize];
                    ja_send[i][j] = jcol;
                }
            }
        }
    }


    numsend = handler->ptrvsend[handler->npsend];
    numrecv = handler->ptrvrecv[handler->nprecv];

    nodv = handler->nodv;
    n = lmat->n + nodv;
    //n = lmat->nc;
    printf(" handler->nodv = %d\n",nodv);

    new_mat->n = n;//new mat, may we need to initialise more members here
    new_mat->bszc = NULL;
    new_mat->D = NULL;
    new_mat->lsize = 0;

    PARMS_NEWARRAY(new_mat->bsz, n+1);
    PARMS_MEMCPY(new_mat->bsz,bszc , n+1);

    PARMS_NEWARRAY(new_mat->nzcount, n);
    PARMS_NEWARRAY(new_mat->ja,     n);
    PARMS_NEWARRAY(new_mat->ba,     n);

    /* copy the local matrix to new_mat */
    for (i = 0; i < lmat->n; i++) {
        dimR = B_DIM(bsz, i);
        nnz = new_mat->nzcount[i] = lmat->nzcount[i];
        if (nnz) {
            PARMS_NEWARRAY(new_mat->ja[i], nnz);
            PARMS_NEWARRAY(new_mat->ba[i], nnz);
            PARMS_MEMCPY(new_mat->ja[i], lmat->ja[i], nnz);
            //PARMS_MEMCPY(new_mat->ba[i], lmat->ba[i], nnz);
            //we may need to use copyBData here
            for (j = 0; j < nnz; ++j)	{
                jcol = lmat->ja[i][j];
                dimC = B_DIM(bszc, jcol);
                PARMS_NEWARRAY(new_mat->ba[i][j], dimC*dimR);
                copyBData(dimC, dimR, new_mat->ba[i][j], lmat->ba[i][j], 0);
            }
        }
    }

    if (handler->npsend+handler->nprecv) {
        PARMS_NEWARRAY(sstatus, handler->npsend+handler->nprecv);
        PARMS_NEWARRAY(sreq,    handler->npsend+handler->nprecv);
        rstatus = sstatus + handler->npsend;
        rreq    = sreq + handler->npsend;
    }

    tag = 800;

    /* exchange external matrix */
    if (numrecv) {
        PARMS_NEWARRAY(rbuf, numrecv);
    }

    /* receive the number of entries for each row in the extended
     matrix */
    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        MPI_Irecv(&rbuf[pos], length, MPI_INT, handler->procs_recv[i],
                  tag, comm, &rreq[i]);//int index exchange
    }

    if (numsend) {//remain the same
        PARMS_NEWARRAY(sbuf, numsend);//sbuf is to store the nnz of the rows
    }


    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            sbuf[pos+j] = ia_send[lindex];
        }
        MPI_Isend(&sbuf[pos], length, MPI_INT, handler->procs_send[i],
                  tag, comm, &sreq[i]);
    }
    MPI_Waitall(handler->nprecv, rreq, rstatus);
    MPI_Waitall(handler->npsend, sreq, sstatus);

    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        for (j = 0; j < length; j++) {
            new_mat->nzcount[lmat->n+pos+j] = rbuf[pos+j];
        }
    }

    j = 0;
    for (i = 0; i < nodv; i++) {
        length = new_mat->nzcount[i+lsize];
        j += length;
        if (length) {
            PARMS_NEWARRAY(new_mat->ja[i+lmat->n], length);
            PARMS_NEWARRAY(new_mat->ba[i+lmat->n], length);
        }
    }
    if (j) {// j is the sum of length of each row in the receiving buffer
        PARMS_NEWARRAY(ja_rbuf, j);//int index
        PARMS_NEWARRAY(ba_rbuf, j);//FLOAT value, a array for all the rows
        PARMS_NEWARRAY(ba_rbuf_nB, j);//new
    }

    PARMS_NEWARRAY(pia_send, handler->npsend+1);
    PARMS_NEWARRAY(pia_recv, handler->nprecv+1);
    pia_send[0] = 0;
    pia_recv[0] = 0;


    /* receive column indices of each row */
    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            m = new_mat->nzcount[lmat->n+pos+j];
            k += m;
        }
        pia_recv[i+1] = pia_recv[i] + k;
        MPI_Irecv(&ja_rbuf[pia_recv[i]],  k, MPI_INT,
                handler->procs_recv[i], tag, comm, &rreq[i]);
    }

    k = 0;
    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            k += m;
        }
    }

    if (k) {// k is the sum of length of each row in the sending buffer
        PARMS_NEWARRAY(ja_sbuf, k);//int index
        PARMS_NEWARRAY(ba_sbuf, k);//FLOAT value, a array for all the rows
        PARMS_NEWARRAY(ba_sbuf_nB, k);//new
    }

    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            PARMS_MEMCPY(&ja_sbuf[pia_send[i]+k], ja_send[lindex], m);
            k += m;
        }
        pia_send[i+1] = pia_send[i] + k;//pia_send is like ptrvsend, for packed all row's data
        MPI_Isend(&ja_sbuf[pia_send[i]], k, MPI_INT, handler->procs_send[i],
                tag, comm, &sreq[i]);
    }
    MPI_Waitall(handler->nprecv, rreq, rstatus);
    MPI_Waitall(handler->npsend, sreq, sstatus);

    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            m = new_mat->nzcount[lmat->n+pos+j];
            PARMS_MEMCPY(new_mat->ja[lmat->n+pos+j], &ja_rbuf[pia_recv[i]+k],
                    m);
            k += m;
        }
    }

    /*construct pia_send_b and pia_recv_b, they are head index for unpacked data*/
    /*in the traversal, we also inserted value into ba_rbuf_nB and ba_sbuf_nB*/
    PARMS_NEWARRAY(pia_send_b, handler->npsend+1);
    PARMS_NEWARRAY(pia_recv_b, handler->nprecv+1);
    pia_send_b[0] = 0;
    pia_recv_b[0] = 0;


    /*for receiving*/
    knz = 0;//calculate the head index for ba_rbuf_nB
    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            index = lmat->n+pos+j;
            nnz = new_mat->nzcount[index];
            dimR = B_DIM(bszc, index);
            if(nnz){
                for (jj = 0; jj < nnz; jj++){
                    jcol = new_mat->ja[index][jj];//here jcol is the global index
                    dimC = is->nB[jcol];//so use is->nB[jcol] to get the block size
                    blocksz = dimR*dimC;
                    k += blocksz;
                    indexnB = knz + jj;
                    ba_rbuf_nB[indexnB] = blocksz;
                }
            }
            knz += nnz;
            //k += m;
        }
        pia_recv_b[i+1] = pia_recv_b[i] + k;
        //MPI_Irecv(&ja_rbuf[pia_recv[i]],  k, MPI_INT,
        //      handler->procs_recv[i], tag, comm, &rreq[i]);
    }


    if (pia_recv_b[handler->nprecv]) {
        //PARMS_NEWARRAY(ja_sbuf, pia_recv_b[i]);
        PARMS_NEWARRAY(ba_rbuf_b, pia_recv_b[handler->nprecv]);
    }

    /*for sending */
    knz = 0;
    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            index = handler->vlist_send[pos+j];
            nnz = lmat->nzcount[index];
            dimR = B_DIM(bsz, index);

            if(nnz){
                for (jj = 0; jj < nnz; jj++){
                    jcol = lmat->ja[index][jj];
                    dimC = B_DIM(bszc, jcol);
                    blocksz = dimR*dimC;
                    k += blocksz;
                    indexnB = knz + jj;
                    ba_sbuf_nB[indexnB] = blocksz;
                }
            }
            knz += nnz;
            /* pindex = parms_TableGet(is->vstable, lindex); */
            /* lindex = *pindex; */
            /* m = ia_send[lindex]; */
            //PARMS_MEMCPY(&pj_sbuf[pia_send[i]+k], ja_send[lindex], m);
            //k += m;
        }
        pia_send_b[i+1] = pia_send_b[i] + k;//pia_send is like ptrvsend, for packed all row's data
        /* MPI_Isend(&pj_sbuf[pia_send[i]], k, MPI_INT, handler->procs_send[i],  */
        /* 	      tag, comm, &sreq[i]); */
    }

    if (pia_send_b[handler->npsend]) {
        //PARMS_NEWARRAY(ja_sbuf, pia_recv_b[i]);
        PARMS_NEWARRAY(ba_sbuf_b, pia_send_b[handler->npsend]);
    }

    /* receive values of each row */


#if defined(DBL_CMPLX)  
    for (i = 0; i < handler->nprecv; i++) {
        k = pia_recv_b[i+1] - pia_recv_b[i];//use pia_recv_b instead of pia_recv
        MPI_Irecv(&ba_rbuf_b[pia_recv_b[i]], k, MPI_CMPLX,
                handler->procs_recv[i], tag, comm, &rreq[i]);
    }

    nBsum = 0;//to get the head index of ba_sbuf_b to store the unpacked data
    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            PARMS_MEMCPY(&ba_sbuf[pia_send[i]+k], a_send[lindex], m);
            k += m;
        }

        kk = pia_send_b[i+1] - pia_send_b[i];//length of each sending
        for (j = pia_send[i]; j < pia_send[i+1]; j++) {//unpack data
            PARMS_MEMCPY(&ba_sbuf_b[nBsum], ba_sbuf[j], ba_sbuf_nB[j]);
            nBsum += ba_sbuf_nB[j];
        }
        //MPI_Isend(&ba_sbuf_b[pia_send_b[i]], k, MPI_CMPLX,
        //      handler->procs_send[i], tag, comm, &sreq[i]);
        MPI_Isend(&ba_sbuf_b[pia_send_b[i]], kk, MPI_CMPLX,
                handler->procs_send[i], tag, comm, &sreq[i]);
    }
#else 
    for (i = 0; i < handler->nprecv; i++) {
        k = pia_recv_b[i+1] - pia_recv_b[i];//use pia_recv_b instead of pia_recv
        MPI_Irecv(&ba_rbuf_b[pia_recv_b[i]], k, MPI_DOUBLE,
                handler->procs_recv[i], tag, comm, &rreq[i]);
    }

    nBsum = 0;//to get the head index of ba_sbuf_b to store the unpacked data
    for (i = 0; i < handler->npsend; i++) {
        pos = handler->ptrvsend[i];
        length = handler->ptrvsend[i+1]-handler->ptrvsend[i];
        k = 0;
        for (j = 0; j < length; j++) {
            lindex = handler->vlist_send[pos+j];
            pindex = parms_TableGet(is->vstable, lindex);
            lindex = *pindex;
            m = ia_send[lindex];
            PARMS_MEMCPY(&ba_sbuf[pia_send[i]+k], a_send[lindex], m);//???
            k += m;
        }

        kk = pia_send_b[i+1] - pia_send_b[i];//length of each sending
        for (j = pia_send[i]; j < pia_send[i+1]; j++) {//unpack data
            PARMS_MEMCPY(&ba_sbuf_b[nBsum], ba_sbuf[j], ba_sbuf_nB[j]);
            nBsum += ba_sbuf_nB[j];
        }

        MPI_Isend(&ba_sbuf_b[pia_send_b[i]], kk, MPI_DOUBLE,
                handler->procs_send[i], tag, comm, &sreq[i]);
    }
#endif

    MPI_Waitall(handler->nprecv, rreq, rstatus);
    MPI_Waitall(handler->npsend, sreq, sstatus);

    MPI_Comm_free(&comm);

    /*pack ba_rbuf_b data (double) to ba_rbuf (BData) type*/

    nBsum = 0;
    for (i = 0; i < pia_recv[handler->nprecv]; i++){
        ba_rbuf[i] = &ba_rbuf_b[nBsum];
        nBsum += ba_rbuf_nB[i];
    }

    //output_intvectorpa("ba_rbuf_nB",ba_rbuf_nB,0, pia_recv[handler->nprecv]);
    //output_dblvectorpa("ba_rbuf_b",ba_rbuf_b,0, pia_recv_b[handler->nprecv]);


    /* for (i = 0; i < pia_recv[handler->nprecv]; i++){ */
    /*   PARMS_NEWARRAY(ba_rbuf[i], ba_rbuf_nB[i]); */
    /*   PARMS_MEMCPY(ba_rbuf[i], &ba_rbuf_b[nBsum], ba_rbuf_nB[i]); */
    /*   //copyBData(ba_rbuf_nB[i], 1, ba_rbuf[i], &ba_rbuf_b[nBsum], 0); */
    /*   //ba_rbuf[i] = &ba_rbuf_b[nBsum]; */
    /*   nBsum += ba_rbuf_nB[i]; */
    /* } */

    for (i = 0; i < handler->nprecv; i++) {
        pos = handler->ptrvrecv[i];
        length = handler->ptrvrecv[i+1] - handler->ptrvrecv[i];
        k = 0;
        for (j = 0; j < length; j++) {
            m = new_mat->nzcount[lmat->n+pos+j];
            PARMS_MEMCPY(new_mat->ba[lmat->n+pos+j], &ba_rbuf[pia_recv[i]+k],
                    m);
            k += m;
        }
    }

    /* free temporary arrays */
    if (handler->npsend) {
        PARMS_FREE(ja_sbuf);
        PARMS_FREE(ba_sbuf);
        PARMS_FREE(ba_sbuf_b);//new
        PARMS_FREE(ba_sbuf_nB);//new
    }



    if (numrecv) {
        PARMS_FREE(rbuf);
    }

    if (numsend) {
        PARMS_FREE(sbuf);
    }


    PARMS_FREE(pia_recv);
    PARMS_FREE(pia_send);
    PARMS_FREE(pia_recv_b);//new
    PARMS_FREE(pia_send_b);//new

    if (handler->npsend+handler->nprecv) {
        PARMS_FREE(sreq);
        PARMS_FREE(sstatus);
    }

    if (is->ninf_send) {
        for (i = 0; i < is->ninf_send; i++) {
            if (ia_send[i]) {
                PARMS_FREE(ja_send[i]);
                PARMS_FREE(a_send[i]);
            }
        }

        PARMS_FREE(ia_send);
        PARMS_FREE(ja_send);
        PARMS_FREE(a_send);
    }

    /* change the global column indices into local ones */
    for (i = lmat->n; i < n; i++) {
        nnz = new_mat->nzcount[i];
        rowj = new_mat->ja[i];
        rowa = new_mat->ba[i];
        k = 0;
        for (j = 0; j < nnz; j++) {
            jcol = rowj[j];
            pindex = parms_TableGet(is->table, jcol);
            if (pindex) {
                if (*pindex < lsize) {
                    rowj[j] = is->perm[*pindex] - start;
                }
                else if (*pindex >= lsize) {
                    rowj[j] = *pindex - start;
                }
                if (k < j) {
                    rowj[k] = rowj[j];
                    rowa[k] = rowa[j];
                }
                k++;
            }
        }
        new_mat->nzcount[i] = k;
        //if (myid == 0)
        //printf("k = %d,  pid = %d \n", k, myid);

        if (k) {// now we have the final row size, so we start to copy the final block data to new_mat->ba[i]
            dimR = B_DIM(bszc, i);
            for (j = 0; j < k; j++){
                jcol = new_mat->ja[i][j];
                dimC = B_DIM(bszc, jcol);
                blocksz = dimR*dimC;
                //	if (myid == 0)
                //printf("k = %d, blocksz = %d, pid = %d \n", k, blocksz, myid);

                PARMS_NEWARRAY(blockbuffer, blocksz);
                PARMS_MEMCPY(blockbuffer, new_mat->ba[i][j], blocksz);
                new_mat->ba[i][j] = blockbuffer;
            }
            PARMS_RESIZE(new_mat->ja[i], k);
            PARMS_RESIZE(new_mat->ba[i], k);
        }
    }

    /* free temporary arrays */

    if (handler->nprecv) {
        PARMS_FREE(ja_rbuf);
        PARMS_FREE(ba_rbuf);
        PARMS_FREE(ba_rbuf_b);//new
        PARMS_FREE(ba_rbuf_nB);//new
    }

    *ext_mat = new_mat;
    *nn = n;
    return 0;
}

static int MatFreeSubMat_b_dvcsr(parms_Mat self, void *mat)
{
    parms_bvcsr vmat;
    int i, nnz, isalloc, j;

    isalloc = self->isalloc;

    //printf("in MatFreeSubMat_b_dvcsr");
    if(isalloc){
        vmat = (parms_bvcsr)mat;
        PARMS_FREE(vmat->bsz);
        for (i = 0; i < vmat->n; i++) {
            nnz = vmat->nzcount[i];
            if (nnz) {
                PARMS_FREE(vmat->ja[i]);
                for( j = 0; j < nnz; j++ ) {
                    PARMS_FREE( vmat->ba[i][j] );
                }
                PARMS_FREE(vmat->ba[i]);
            }
        }

        if (vmat->n) {
            PARMS_FREE(vmat->nzcount);
            PARMS_FREE(vmat->ja);
            PARMS_FREE(vmat->ba);
        }

        PARMS_FREE(vmat);

        mat = NULL;
    }

    return 0;
}



static struct parms_Mat_ops parms_matops_b_dvcsr = {
    MatVec_b_dvcsr,
    MatSetup_b_dvcsr,
    MatSetCommType_dvcsr,
    MatMVPY_b_dvcsr,//MatMVPY_dvcsr,
    MatGetDiag_b_dvcsr,//MatGetDiag_dvcsr,
    MatGetSubMat_b_dvcsr,
    MatExtend_b_dvcsr,
    MatVecSchur_b_dvcsr,
    MatFreeSubMat_b_dvcsr,//MatFreeSubMat_dvcsr,
    MatGetHandler_dvcsr
};



int parms_MatFree_b_dvcsr(parms_Mat *self)
{
    parms_dvcsr   data;
    parms_bvcsr    b_diag_mat, b_offd_mat;

    data = (parms_dvcsr)(*self)->data;

    parms_CommFree(&data->mvhandler);
    b_diag_mat = data->b_diag_mat;
    PARMS_FREE(b_diag_mat->nzcount);
    PARMS_FREE(b_diag_mat->ba);
    PARMS_FREE(b_diag_mat->ja);
    PARMS_FREE(b_diag_mat->bsz);

    //printf("b_diag_mat->D= %p\n",b_diag_mat->D);
    if (b_diag_mat->D) {
        PARMS_FREE(b_diag_mat->D);
    }
    PARMS_FREE(b_diag_mat);

    b_offd_mat = data->b_offd_mat;
    if (b_offd_mat->n) {
        PARMS_FREE(b_offd_mat->nzcount);
        PARMS_FREE(b_offd_mat->ba);
        PARMS_FREE(b_offd_mat->ja);
        PARMS_FREE(b_offd_mat->bsz);
        if (b_offd_mat->D) {
            PARMS_FREE(b_offd_mat->D);
        }
        if (b_offd_mat->bszc) {
            PARMS_FREE(b_offd_mat->bszc);
        }
    }

    //printf("b_offd_mat->bszc = %p\n",b_offd_mat->bszc);

    PARMS_FREE(b_offd_mat);

    PARMS_FREE(data);
    return 0;
}


int parms_MatCreate_b_dvcsr(parms_Mat self)
{
    parms_dvcsr  data;
    parms_Comm   mvhandler;
    MPI_Comm     comm;

    PARMS_MEMCPY(self->ops, &parms_matops_b_dvcsr, 1);
    PARMS_NEW(data);
    comm = self->is->comm;
    parms_CommCreate(&mvhandler, comm);
    data->mvhandler = mvhandler;

    self->data = data;
    return 0;
}
