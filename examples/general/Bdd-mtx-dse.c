#define _GNU_SOURCE
/*---------------------------------------------------------------------- 
 *                           Program dd-HB-dse
 *----------------------------------------------------------------------
 *
 *  In this test program, each processor reads the whole matrix
 *  from file. The matrix is assumed to be in Harwell-Boeing format.
 *  Matrix graph is then  partitioned  using  DSE, a simple partitioning
 *  routine, and scatters the local matrices to each processor. Once
 *  these submatrices are received each processor solves the problem
 *  using preconditioned FGMRES preconditioned with :
 *                       BJ, RAS, SCHUR
 *--------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(__ICC)
#include <mathimf.h>
#else
#include <math.h> 
#endif 
#include "parms.h"
#include "aux.h"
#include "mmio.h"
//#include <time.h>//add

/* typedef struct parms_dvcsr { */
/*   /\*! \var diag_mat The diagonal matrix stored in vcsr format. */
/*    *\/ */
/*   parms_vcsr    diag_mat; */
/*  /\*!  \var offd_mat The off-diagonal matrix stored in vcsr format. */
/*   *\/ */
/*   parms_vcsr    offd_mat; */
/*   /\*! \var mvhandler The parms_Comm object for the matrix-vector product.  */
/*    *\/ */
/*   parms_bvcsr   b_diag_mat;//vbsptr */
/*   parms_bvcsr   b_offd_mat;//vbsptr */
/*   parms_Comm    mvhandler;	 */

/* } *parms_dvcsr; */
#include <dlfcn.h>
#include <stdlib.h>

unsigned long long (*get_mem)() = NULL;

unsigned long long get_memory_usage_()
{
    static int once = 0;
    if (!once)
    {
        get_mem = (unsigned long long (*)())dlsym(RTLD_DEFAULT, "get_memory_usage");
        if (!get_mem)
        {
            get_mem = get_memory_usage_;
            printf("To enable memory profiling you need to LD_PRELOAD the malloc_impl library\n");
            once = 1;
            return 0;
        }
        return get_mem();
    }
    return 0;
}

void print_mem(const char* descr)
{
    int myid;
    double value;
    char unit[3];
    strcpy(unit, "B ");
    value = get_mem();

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if (value > 1.0e3) {value*=1.0e-3; strcpy(unit, "kB");}
    if (value > 1.0e3) {value*=1.0e-3; strcpy(unit, "MB");}
    if (value > 1.0e3) {value*=1.0e-3; strcpy(unit, "GB");}
    if (value > 1.0e3) {value*=1.0e-3; strcpy(unit, "TB");}

    printf("%s on proc %d: %f %s\n", descr, myid, value, unit);
}


int main(int argc, char *argv[])
{

    /* declarations related to Harwell-boeing format for reading the HB
     matri. Second part is related to I/O parameters */
    char mname[MAX_MAT][MAX_LINE], key[MAX_LINE];//, type[3];//guesol[2], title[72],
    int nrhs, nc, n, nnz, mat;//tmp0,tmp2, tmp3, tmp,
    int myid, ierr, i, nloc;
    /* memory usage of the preconditioning matrix */
    double ratio, Bdensity = 0.0;
    /* working array for reading matrix */
    double norm, res1, tpc, ttol;
    FLOAT *a, *rhstmp;
    int *ja, *ia;
    int npro,its, *im, *bim;//bim is for block case
    fprm prm;
    FILE *fp=NULL, *fout=NULL, *mtxfile=NULL;
    char *name, *iluname, *curmat, *currhs, *curname = NULL, buf[40];
    int *nzding = NULL;
    int nBlock, *nB = NULL, *perm = NULL;
    double tib1, tib2, tib3, tib4 = 0.0, blocksize = 0.0;

    /*-------------------- variables related to dse partitioning */
    int *riord, *dom = NULL, *idom = NULL, *mask, *jwk, *link;

    /*-------------------- variables related to pARMS */
    parms_Map       map;
    FLOAT       *x, *y, *rhs, *resvec;
    parms_Mat       A;
    parms_PC        pc;
    parms_Solver    solver;
    parms_Timer     tm;



    /* Viewer object for solver */
    //  parms_Viewer  sv;

    /*-------------------- initialize MPI environment */
    MPI_Init(&argc, &argv);//initilization
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);//return the pid
    MPI_Comm_size(MPI_COMM_WORLD, &npro);//return the number of the processors
    get_mem = get_memory_usage_;


//    tmp0 = 0;
    nrhs = 0;

    /*-------------------- read matrix name from input file */
    prm = malloc(sizeof(*prm));
    if (prm == NULL) {
        fprintf(stderr, "cannot allocate memory for prm\n");
        MPI_Abort(MPI_COMM_WORLD, 66);
    }

    if (argc >= 2) {
        read_param(argv[1],mname, prm);
    }
    else if (argc == 1){
        read_param("inputs",mname, prm);
    }


    /* variable "mname" stores the name of the file in HB format. Read a
     Harwell-Boeing matrix. using wreadmtc c-version of sparsit
     routine - call wreadmtc a first time to determine sizes of
     arrys. read in values on the second call.
  */

    /* --- Begin loop over matrices ----*/
    mat = 0;
    while(mname[mat][1] != '#'){
        a = NULL; ja = NULL; ia = NULL; rhstmp = NULL;
        curmat = strtok(mname[mat], " ");
        if (curmat == NULL)
            continue;
        currhs = strtok(NULL, " ");
        parms_TimerCreate(&tm);
        parms_TimerReset(tm);

        if (myid == 0)
            printf("Reading matrix %s\n", curmat);

        //~ #if defined(DBL_CMPLX)
        //~ zreadmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
        //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
        //~ #else
        //~ readmtc_(&tmp0,&tmp0,&tmp0,mname[mat],a,ja,ia,rhstmp,&nrhs,
        //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
        //~ #endif

        if ((mtxfile = fopen(curmat, "r")) == NULL) {
            fprintf(stderr, "Error opening matrix file\n");
            fprintf(stderr, "filename = %s\n", curmat);
            MPI_Finalize();
            exit(1);
        }

        MM_typecode matcode;
        ierr = mm_read_banner(mtxfile, &matcode);
        if(ierr != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "cannot read banner\n");
            fprintf(stderr, "filename = %s\n", curmat);
            MPI_Finalize();
            exit(1);
        }

        if ((ierr = mm_read_mtx_crd_size(mtxfile, &n, &nc, &nnz, matcode)) != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "Error reading matrix dimensions\n");
            MPI_Finalize();
            exit(1);
        }

        if (myid == 0)
        {
            a = malloc(nnz*sizeof(*a));
            if (a == NULL) {
                fprintf(stderr, "cannot allocate memory for a\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            ja = malloc(nnz*sizeof(*ja));
            if (ja == NULL) {
                fprintf(stderr, "cannot allocate memory for ja\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            ia = malloc(nnz*sizeof(*ia));
            if (ia == NULL) {
                fprintf(stderr, "cannot allocate memory for ia\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
//            tmp2 = n;
//            tmp3 = nnz;

            /*-------------------- Array sizes determined. Now call
                         wreadmtc again for really reading */
            /*
  wreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],&len,a,ja,ia,rhstmp,&nrhs,
        guesol,&n,&nc,&nnz,title,key,type,&ierr);
*/
            //~ #if defined(DBL_CMPLX)
            //~ zreadmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
            //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
            //~ #else
            //~ readmtc_(&tmp2,&tmp3,&tmp,mname[mat],a,ja,ia,rhstmp,&nrhs,
            //~ guesol,&n,&nc,&nnz,title,key,type,&ierr);
            //~ #endif

            ierr = mm_read_mtx_crd_data(mtxfile, n, nc, nnz, ia, ja, a, matcode);
            fclose(mtxfile);

            if ((mtxfile = fopen(curmat, "r")) == NULL) {
                fprintf(stderr, "Error opening matrix file\n");
                fprintf(stderr, "filename = %s\n", curmat);
                MPI_Finalize();
                exit(1);
            }

            MM_typecode matcode;
            ierr = mm_read_banner(mtxfile, &matcode);
            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read banner\n");
                fprintf(stderr, "filename = %s\n", curmat);
                MPI_Finalize();
                exit(1);
            }

            if ((ierr = mm_read_mtx_crd_size(mtxfile, &n, &nc, &nnz, matcode)) != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "Error reading matrix dimensions\n");
                MPI_Finalize();
                exit(1);
            }

            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read matrix\n");
                fprintf(stderr, "filename = %s\n", mname[mat]);
                MPI_Finalize();
                exit(1);
            }

            if(myid == 0){
                if(argc == 3) {
                    if (NULL == (fp = fopen(argv[2], "w"))) {
                        fprintf(stderr, "cannot open file %s\n", argv[2]);
                        exit(1);
                    }
                }
                else {
                    fp = stdout;
                }
                curname = strtok(curmat, "/");
                while (curname != NULL)
                {
                    strcpy(key, curname);
                    curname = strtok(NULL, "/");
                }
                curname = strtok(key, ".");
                fprintf(fp, "\nMatrix %d: %s \n",(mat+1), curname);
                fprintf(fp, "n = %d, nnz = %d\n", n, nnz);
            }

            /* /\*-------------Convert from COO to CSR format ------------*\/ */
            /*     int *jb, *ib; */
            /*     FLOAT *b; */
            /*     b   = malloc(nnz*sizeof(*b)); */
            /*     if (b == NULL) { */
            /*       fprintf(stderr, "cannot allocate memory for b\n"); */
            /*       MPI_Abort(MPI_COMM_WORLD, 66); */
            /*     } */
            /*     jb  = malloc(nnz*sizeof(*jb)); */
            /*     if (jb == NULL) { */
            /*       fprintf(stderr, "cannot allocate memory for jb\n"); */
            /*       MPI_Abort(MPI_COMM_WORLD, 66); */
            /*     } */
            /*     ib  = malloc((n+1)*sizeof(*ib)); */
            /*     if (ib == NULL) { */
            /*       fprintf(stderr, "cannot allocate memory for ib\n"); */
            /*       MPI_Abort(MPI_COMM_WORLD, 66); */
            /*     } */
            /*     job = 1; */


            output_intvector("ia.coo",ia,0, nnz);
            output_intvector("ja.coo",ja,0, nnz);
            output_dblvector("a.coo",a,0, nnz);



            /* #if defined(DBL_CMPLX) */
            /*     zcoocsr_(&n, &nnz, a, ia, ja, b, jb, ib); */
            /* #else */
            /*     coocsr_(&n, &nnz, a, ia, ja, b, jb, ib); */
            /* #endif */
            /*---------------Free COO matrix ------------------*/
            /* /\*-------------Convert from COO to C-CSR format ------------*\/ */

            csptr csmat = NULL;

            //mat = malloc(sizeof(*(csptr)));
            csmat = malloc(sizeof(*csmat));

            coo2csptr(n, nnz, a, ia, ja, csmat);
            /* outputcsmat(mat,"mat.coo",1);  */

            /* exit(1); */

            free(ia);
            free(a);
            free(ja);

            /* ia = ib; */
            /* ja = jb; */
            /* a = b; */

            //   double t1;
            //   t1 = sys_timer();
            //csptr csmat = NULL;
            //ddouble aa;
            vbsptr vbmat = NULL;
            //BData bb;
            //p4ptr prev;
            /* csmat = (csptr)Malloc( sizeof(SparMat), "testdriver" );//csmat = (csptr)malloc(sizeof(SparMat)); */
            /* colunms2csptr(n, ia, ja, a, csmat); //int colunms2csptr(int n, int *ia, int *ja, FLOAT *a, csptr mat)  */
            //if (myid == 0)
            /* outputcsmatpa(csmat,"csmat.coo",1);  */
            /* exit(1); */
            /*--------------------create a timer */
            //~
            //~ Cliques *MC = cliques_new();
            //~ Clique *R = clique_new();
            //~ Clique *P = clique_new();
            //~ for (i = 0; i < n; ++i)
            //~ clique_append(P, i);
            //~ Clique *X = clique_new();
            //~ get_maximal_cliques(MC, csmat, R, P, X);
            //~
            //~ clique_free(R);
            //~ clique_free(P);
            //~ clique_free(X);
            //~ cliques_free(MC);
            //~
            //~ MPI_Abort(MPI_COMM_WORLD, 66);

            printf("prm->eps = %f\n",prm->eps);
            //~ ierr = init_blocks_cliques( csmat, &nBlock, &nB, &perm, prm->eps);
            if (prm->cosine)
                ierr = init_blocks( csmat, &nBlock, &nB, &perm, prm->eps);//int init_blocks( csptr csmat, int *pnBlock, int **pnB, int **pperm, double eps)//parms_PCSetup(pc);
            else
                ierr = init_blocks_density( csmat, &nBlock, &nB, &perm, prm->eps);
            //printf("nBlock: %d\n", nBlock);getchar();
            tib1 =  parms_TimerGet(tm);
            printf("\ntime on init=%f\n",tib1);
            if(ierr != 0) {
                //fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "*** in init_blocks ierr != 0 ***\n");
                MPI_Finalize();
                exit(1);
            }
            //output_intvector("perm.coo",perm,0, n);getchar();
            if( dpermC( csmat, perm ) != 0 ) {
                fprintf( stderr, "*** dpermC error ***\n" );
                MPI_Finalize();
                exit(1);
            }
            tib2 =  parms_TimerGet(tm);
            printf("\ntime on dpermC=%f\n",tib2-tib1);
            /*-------------------- convert to block matrix. */
            vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
            //ierr = csrvbsrC( 1, nBlock, nB, csmat, vbmat );
            ierr = csrvbsrC_new( 1, nBlock, nB, csmat, vbmat );
            //outputvbmatpa(vbmat,"vbmat2.coo",1);//int outputvbmat( vbsptr vbmat, char *filename, int onebase)
            tib3 =  parms_TimerGet(tm);
            printf("\ntime on csrvbsrC_new=%f\n",tib3-tib2);
            //~ outputvbmatpa(vbmat,"vbmat2.coo",1);


            //exit(1);
            //if (myid == 0)
            //outputvbmatpa(vbmat,"vbmat.coo",1);//int outputvbmat( vbsptr vbmat, char *filename, int onebase)

            //sprintf(buf, "%8.2fs", tib);
            //sprintf(buf, "%8.2fs", tpc);
            //fprintf(fp, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
            //fprintf(fp, "The time cost for init_blocks %7s %-s\n", "=", strtok(buf, " "));
            blocksize = (double)csmat->n / (double)nBlock;
            Bdensity = (double)nnzCS( csmat ) / (double)memVBMat( vbmat ) * 100;
            if (myid == 0)
            {
                printf("\n Bsize=%-7f, Bdensity=%-7f\n",blocksize, Bdensity);
                //~ FILE* bfile = fopen("bsize.txt", "w");
                //~ fprintf(bfile, "%f\n", Bdensity);
                //~ fclose(bfile);
            }
            //~ MPI_Abort(MPI_COMM_WORLD, 0);
            //output_dblvectorpa("rhstmpp",rhstmpp, 0, n);

            //cleanCS(csmat);//int cleanCS(csptr amat)
            int nbb, *bia, *bja;
            BData *ba;
            bia = (int*)malloc((vbmat->n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
            bja = (int*)malloc(nnz*sizeof(int));
            ba = (BData*)malloc(nnz*sizeof(BData));
            vbsptr2colunms(vbmat, &nbb, bia, bja, ba);//int vbsptr2colunms(vbsptr mat, int *n, int *ia, int *ja)
            //iia = (int*)malloc(n*sizeof(int)   printf("\n Bsize=%-7f\n",blocksize););//w =(double*)malloc(n*sizeof(double));
            //jja = (int*)malloc(nnz*sizeof(int));
            //csptr2colunms(csmat, &nn, iia, jja);//int csptr2colunms(csptr mat, int *n, int *ia, int *ja)//
            bja = (int*)realloc(bja, nbb*sizeof(int));
            //printf("bja[0]=%d",bja[0]);
            //output_intvector("bia.coo",bia,0, nBlock);getchar();
            //output_intvector("bja.coo",bja,0, nbb);getchar();


            tib4 =  parms_TimerGet(tm);
            printf("\ntime on csrvbsrC_new=%f, time on total process = %f\n",tib4-tib3, tib4);
            /*--------------------get the elapsed time spent on creating PC */

            //int dpermC(csptr mat, int *perm)


            //printf("\n npro=%d\n",npro);
            //int nn, *iia, *jja;
            //iia = (int*)malloc((n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
            //jja = (int*)malloc(nnz*sizeof(int));
            //csptr2colunms(csmat, &nn, iia, jja);//int csptr2colunms(csptr mat, int *n, int *ia, int *ja)//
            //  im = (int *)malloc(n*sizeof(int));
            //printf("nn = %d\n", nn );getchar();
            //output_intvector("iia.coo",iia,0, nn+1);getchar();
            //output_intvector("jja.coo",jja,0, nnz);getchar();




            idom = malloc((npro+1)*sizeof(*idom));
            if (idom == NULL) {
                fprintf(stderr, "cannot allocate memory for idom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            dom = malloc(nBlock*sizeof(*dom));
            if (dom == NULL) {
                fprintf(stderr, "cannot allocate memory for dom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }

            if (npro == 1) {
                for (i = 0; i < nBlock; i++) {
                    dom[i] = i+1;
                }
                idom[0] = 1;
                idom[1] = nBlock+1;
            }
            else {
                riord = malloc(nBlock*sizeof(*riord));
                if (riord == NULL) {
                    fprintf(stderr, "cannot allocate memory for riord\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                mask = malloc(nBlock*sizeof(*mask));
                if (mask == NULL) {
                    fprintf(stderr, "cannot allocate memory for mask\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                jwk = malloc(2*nBlock*sizeof(*jwk));
                if (jwk == NULL) {
                    fprintf(stderr, "cannot allocate memory for jwk\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                link = malloc(nBlock*sizeof(*link));
                if (link == NULL) {
                    fprintf(stderr, "cannot allocate memory for link\n");
                    MPI_Abort(MPI_COMM_WORLD, 66);
                }
                dse_(&nBlock, bja, bia, &npro, riord, dom, idom, mask, jwk, link);
                free(riord);
                free(mask);
                free(jwk);
                free(link);
            }

            free(ba);
            free(bja);
            free(bia);

            nzding = (int *) calloc(npro, sizeof(int));
            int i1, i3;
            int i2 = 0, i4;
            for (i1 = 0; i1 < nBlock; ++i1)
            {
                if (i1 == idom[i2+1]-1)
                    i2++;
                i4 = dom[i1]-1;
                for (i3 = vbmat->bsz[i4]; i3 < vbmat->bsz[i4+1]; ++i3)
                    nzding[i2] += csmat->nnzrow[i3];
            }

            print_mem("end of proc 1");
            cleanCS( csmat );
            cleanVBMat( vbmat );//int cleanVBMat( vbsptr vbmat )
        }
        tib2 =  parms_TimerGet(tm);
        MPI_Bcast(&nBlock, 1, MPI_INT, 0, MPI_COMM_WORLD);
        if (myid != 0)
        {
            nB = (int *) malloc(nBlock * sizeof(int));
            if (nB == NULL) {
                fprintf(stderr, "cannot allocate memory for nB\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            perm = (int *) malloc(n * sizeof(int));
            if (perm == NULL) {
                fprintf(stderr, "cannot allocate memory for perm\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            idom = malloc((npro+1)*sizeof(*idom));
            if (idom == NULL) {
                fprintf(stderr, "cannot allocate memory for idom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            dom = malloc(nBlock*sizeof(*dom));
            if (dom == NULL) {
                fprintf(stderr, "cannot allocate memory for dom\n");
                MPI_Abort(MPI_COMM_WORLD, 66);
            }
            nzding = (int *) malloc(npro * sizeof(int));
        }
        MPI_Bcast(nB, nBlock, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(perm, n, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(idom, npro+1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(dom, nBlock, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(nzding, npro, MPI_INT, 0, MPI_COMM_WORLD);
        printf("nzding %d on %d\n", nzding[myid], myid);
        print_mem("after free");

//        output_intvectorpa("idom.coo", idom, 0, npro+1);
//        output_intvectorpa("dom.coo", dom, 0, nBlock);

        tib1 =  parms_TimerGet(tm);
        printf("\ntime on send=%f\n",tib1 - tib2);
        /*
    for (i = 0; i < n; i++) {
      dom[i] = i+1;
    }
    for(i=1; i<npro; i++)
    idom[i] = idom[i-1]+(n/npro);
    idom[npro+1] = n+1;
printf("idom[%d] = %d; dom[%d] = %d \n", myid, idom[myid], idom[myid]-1,dom[idom[myid]-1]);
*/
        /*-------------------- Create map object */
        parms_MapCreateFromPtr(&map, n, dom, idom, MPI_COMM_WORLD, 1, NONINTERLACED);
        //map->nB = nB;
        parms_Map_Assign_blockstructure(map, nB);//int parms_Map_Assign_blockstructure(parms_Map self, int *nB)//
        //output_intvector("map->nB.coo",map->nB,0, nBlock);getchar();


        nloc = parms_MapGetLocalSize(map);

        //output_intvectorpa("mapperm_p",map->lvars,0, nloc);getchar();

        //getlocalsize
        /* Free dom and idom */

        a = malloc(nzding[myid]*sizeof(*a));
        if (a == NULL) {
            fprintf(stderr, "cannot allocate memory for a\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        ja = malloc(nzding[myid]*sizeof(*ja));
        if (ja == NULL) {
            fprintf(stderr, "cannot allocate memory for ja\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
        ia = malloc(nzding[myid]*sizeof(*ia));
        if (ia == NULL) {
            fprintf(stderr, "cannot allocate memory for ia\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }

//        ierr = mm_partial_read_mtx_crd_data(mtxfile, n, nc, nnz, ia, ja, a, idom, dom, perm, nB, nBlock, matcode);
        ierr = mm_partial_read_mtx_crd_data_new(mtxfile, n, nc, nnz, ia, ja, a, idom, dom, perm, nB, nBlock, matcode);
        fclose(mtxfile);

        printf("n value is %d, nc value is %d\n", n, nc);//%f %p %s %c


        free(dom);
        free(idom);

        if(ierr != 0) {
            fprintf(stderr, "ierr = %d\n", ierr);
            fprintf(stderr, "cannot read matrix\n");
            fprintf(stderr, "filename = %s\n", mname[mat]);
            MPI_Finalize();
            exit(1);
        }

        rhstmp = malloc(n*sizeof(*rhstmp));
        if (rhstmp == NULL) {
            fprintf(stderr, "cannot allocate memory for rhstmp\n");
            MPI_Abort(MPI_COMM_WORLD, 66);
        }
//        if(nrhs != 0)
//            tmp = 3;
//        else
//            tmp = 2;

        if (currhs != NULL)
        {
            if ((mtxfile = fopen(currhs, "r")) == NULL) {
                fprintf(stderr, "Error opening rhs file\n");
                MPI_Finalize();
                exit(1);
            }

            ierr = mm_read_banner(mtxfile, &matcode);
            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read banner\n");
                fprintf(stderr, "filename = %s\n", currhs);
                MPI_Finalize();
                exit(1);
            }

            ierr = mm_read_mtx_array_size(mtxfile, &n, &nc);

            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read array size\n");
                fprintf(stderr, "filename = %s\n", currhs);
                MPI_Finalize();
                exit(1);
            }

            nrhs = 1;
            ierr = mm_read_array_data(mtxfile, n, rhstmp, matcode);
            fclose(mtxfile);

            if(ierr != 0) {
                fprintf(stderr, "ierr = %d\n", ierr);
                fprintf(stderr, "cannot read rhs\n");
                fprintf(stderr, "filename = %s\n", currhs);
                MPI_Finalize();
                exit(1);
            }
        }
        else nrhs = 0;

        csptr csmat = NULL;

        //mat = malloc(sizeof(*(csptr)));
        csmat = malloc(sizeof(*csmat));

        coo2csptr(n, nzding[myid], a, ia, ja, csmat);//nzding is the local length array
        printf("nzding[myid] value is %d, myid = %d\n", nzding[myid], myid);//%f %p %s %c

        /* outputcsmat(mat,"mat.coo",1);  */

        /* exit(1); */

        free(ia);
        free(a);
        free(ja);
        vbsptr vbmat = NULL;
        /*--------------------create a timer */
        tib1 =  parms_TimerGet(tm);
        //output_intvector("perm.coo",perm,0, n);getchar();
        if( dpermC( csmat, perm ) != 0 ) {
            fprintf( stderr, "*** dpermC error ***\n" );
            MPI_Finalize();
            exit(1);
        }
        tib2 =  parms_TimerGet(tm);
        printf("\ntime on dpermC=%f\n",tib2-tib1);
        /*-------------------- convert to block matrix. */
        vbmat = (vbsptr)Malloc( sizeof(VBSparMat), "main" );
        //ierr = csrvbsrC( 1, nBlock, nB, csmat, vbmat );
        ierr = csrvbsrC_new( 1, nBlock, nB, csmat, vbmat );
//        outputvbmatpa(vbmat,"vbmatlocal.coo",1);//int outputvbmat( vbsptr vbmat, char *filename, int onebase)
        tib3 =  parms_TimerGet(tm);
        printf("\ntime on csrvbsrC_new=%f\n",tib3-tib2);

        FLOAT *rhstmpp;
        rhstmpp = (FLOAT*)malloc(n*sizeof(FLOAT));

        //output_intvectorpa("perm_init.coo",perm,0, n);
        // output_dblvectorpa("rhstmp",rhstmp, 0, n);
        printf("n value is %d\n", n);//%f %p %s %c


        for( i = 0; i < n; i++ )
            rhstmpp[perm[i]] = rhstmp[i];
        //output_dblvectorpa("rhstmpp",rhstmpp, 0, n);

        //cleanCS(csmat);//int cleanCS(csptr amat)
        int nbb, *bia, *bja;
        BData *ba;
        bia = (int*)malloc((vbmat->n+1)*sizeof(int));//w =(double*)malloc(n*sizeof(double));
        bja = (int*)malloc(nzding[myid]*sizeof(int));
        ba = (BData*)malloc(nzding[myid]*sizeof(BData));
        vbsptr2colunms(vbmat, &nbb, bia, bja, ba);//int vbsptr2colunms(vbsptr mat, int *n, int *ia, int *ja)
        //iia = (int*)malloc(n*sizeof(int)   printf("\n Bsize=%-7f\n",blocksize););//w =(double*)malloc(n*sizeof(double));
        //jja = (int*)malloc(nnz*sizeof(int));
        //csptr2colunms(csmat, &nn, iia, jja);//int csptr2colunms(csptr mat, int *n, int *ia, int *ja)//
        bja = (int*)realloc(bja, nbb*sizeof(int));
        print_mem("end of partial load");

        free(nzding);

        /*-------------------- create a distributed matrix based on map */
        parms_MatCreate(&A, map);
        /*-------------------- Insert values into A */
        /* Initialize global indices array */
        im = (int *)malloc(n*sizeof(int));
        for (i = 0; i < n; i++) {
            im[i] = i+1;
        }

        bim = (int *)malloc(nBlock*sizeof(int));
        for (i = 0; i < nBlock; i++) {
            bim[i] = i+1;
        }
        parms_MatSetValues_b(A, nBlock, bim, bia, bja, ba, INSERT);//parms_MatSetValues(A, n, im, ia, ja, a, INSERT);

        //parms_MatSetValues_b(A, n, im, ia, ja, a, INSERT);//int parms_MatSetValues_b(parms_Mat self, int m, int *im, int *ia, int *ja, BData *values, INSERTMODE mode)
        //  parms_MatSetElementMatrix(A, n, im, ia, ja, a, INSERT);
        //  parms_MatAssembleElementMatrix(A);
        /*-------------------- free the matrix stored in CSR format(a,ja,ia) */
        /* free(a); */
        /* free(ja); */
        /* free(ia); */


        free(ba);
        free(bja);
        free(bia);

        int llsize;
        //printf("llsizetest=%d\n", llsize);
        llsize = parms_MapGetLocallSize(map);
        printf("llsize=%d, nloc = %d\n", llsize, nloc);
        /*-------------------- Create distributed vectors based on map */
        x = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        rhs = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        y = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        resvec = (FLOAT *)malloc(llsize*sizeof(FLOAT));
        //outputcsmatpa(A->aux_data,"aux_data",1);//int outputcsmat ( csptr mat, char *filename, int onebase){
        //output_intvectorpa("localpermbefore",A->is->perm,0, nloc);getchar();


        /*-------------------- Setup the matrix and communication structure */
        parms_Mat_B_Setup(A);//parms_MatSetup(A);
        printf("\n npro=%d\n",npro);

        //MPI_Barrier(MPI_COMM_WORLD);
        //exit(1);
        /*-------------------- Copy rhs or Setup artifical right-hand-side vector */
        //nrhs = 0;
        //output_intvectorpa("localpermafter",A->is->perm,0, nloc);getchar();
        /*  parms_Viewer v; */
        /*   /\* parms_ViewerCreate(&v, "foomatcoo"); *\/ */
        /*   /\* parms_MatView(A, v);//int parms_MatView_vcsr(parms_Mat self, parms_Viewer v) *\/ */
        /*   /\* parms_MatViewCOO(A, v);//int parms_MatViewCOO_vcsr(parms_Mat self, parms_Viewer v) *\/ */
        /*   /\* parms_ViewerFree(&v); *\/ */
        /*   /\* MPI_Barrier(MPI_COMM_WORLD); *\/ */
        /*   /\* exit(1); *\/ */
        /* //to output the diag and off mat of each processor */
        /*   parms_dvcsr data; */
        /*   data    = (parms_dvcsr)A->data; */

        /*   parms_ViewerCreate(&v, "foomvhandler"); */
        /*   parms_CommView(data->mvhandler, v);//int parms_CommView(parms_Comm self, parms_Viewer v) */
        /*   //parms_ViewerFree(&v); */



        /*   parms_ViewerCreate(&v, "foomapaftersetup"); */
        /*   parms_MapView(map, v); */
        /*   //parms_ViewerFree(&v); */


        /*----- permutes right hand side -------------------------------------*/
        //    for( i = 0; i < io.ndim; i++ )
        //      rhs[perm[i]] = rhs0[i];
        printf("nrhs value is %d.\n", nrhs);//%f %p %s %c

        if(!nrhs)
        {
            for(i=0; i<llsize; i++)
            {
                x[i] = 1.0;
            }
            parms_MatVec(A, x, rhs);
        }
        else
        {
            parms_VecSetValues_b(rhs, nBlock, bim, rhstmpp, INSERT, map);//changed//parms_VecSetValues_b(rhs, n, im, rhstmp, INSERT, map);
//            printf(" in arificial vector loading.\n");//%f %p %s %c

        }

        free(rhstmp);
        free(rhstmpp);
        free(im);

        free(bim);

        //output_dblvectorpa("rhs",rhs,0, llsize);getchar();

        /*--------------------Setup initial guess to a vector of all-0. */
        for(i=0; i<llsize; i++)
        {
            x[i] = 0.0;
        }

//                output_dblvector("xvector",x,0, llsize);

        /*--------------------Get 2-norm of initial residual */

        parms_MatVec(A, x, resvec);
//                output_dblvector("resvecdbl.coo",resvec,0, llsize);getchar();
//                output_dblvector("rhs.coo",rhs,0, llsize);getchar();

//        printf("llsize value is %d.\n", llsize);getchar();//%f %p %s %c

        for(i=0; i<llsize; i++)
        {
            resvec[i] = resvec[i] - rhs[i];
//            printf("resvec[i] value is %20.16e.\n", resvec[i]);//%f %p %s %c

        }

//        output_dblvector("resvecdbl.coo",resvec,0, llsize);getchar();
        parms_VecGetNorm2_b(resvec, &norm, map);

        //printf("norm = %20.16e\n",norm);
        free(resvec);
        /*--------------------Create preconditioner based on the matrix A. */
        parms_PCCreate_b(&pc, A);

        /*--------------------set parameters for pc */
        set_pc_params_b(pc, prm);


        /*--------------------reset the timer */
        parms_TimerReset(tm);
        parms_PCSetup_b(pc);


        /*   parms_Viewer v; */
        /*   parms_ViewerCreate(&v, "PC_Bdd");  */
        /*   parms_PCView(pc, v); */

        /*     FLOAT *xx, *yy; */
        /*     PARMS_NEWARRAY(xx, llsize); */
        /*     PARMS_NEWARRAY(yy, llsize); */
        /*     for (i=0; i < llsize; ++i) */
        /*         xx[i] = 1.0; */

        /*     /\* for(i=0;i<llsize;i++) {  *\/ */
        /*     /\*   xx[i]=(double)(rand()/(1000000));  *\/ */
        /*     /\*   //printf("xx[i] = %d \n",j);  *\/ */
        /*     /\* }  *\/ */
        /*     output_dblvectorpa("pcinitial_Bdd",xx, 0, llsize); */

        /*     parms_MatVec(A, xx, yy); */

        /*     //yy = rhs; */

        /*     for (i=0; i < llsize; ++i) */
        /*         xx[i] = 0.0; */
        /*     parms_PCApply(pc, yy, xx); */
        /*     output_dblvectorpa("pcsol_Bdd",xx, 0, llsize); */
        /*     output_dblvectorpa("pcmv_Bdd",yy, 0, llsize); */

//        printf("test fseek\n");
//        MPI_Barrier(MPI_COMM_WORLD);
//        MPI_Finalize();
//        exit(1);


        /*--------------------get the elapsed time spent on creating PC */
        tpc =  parms_TimerGet(tm);

        printf("The time for pc setup %f in proc %d\n", tpc, myid);


        /*--------------------get the ratio of the number of nonzero entries in the
    preconditioning matrix to that in the original matrix */
        parms_PCGetRatio(pc, &ratio);//block version is not ready
        /*--------------------pause the timer */
        parms_TimerPause(tm);



        /*--------------------Create a solver based on A and pc */

        parms_SolverCreate(&solver, A, pc);
        /*--------------------Set the solver type */
        parms_SolverSetType_b(solver, SOLFGMRES);



        /*--------------------set parameters for solver */
        set_solver_params(solver, prm);



        /*--------------------set up solver -- no longer needed - DOK */
        //parms_SolverSetup(solver);

        /*--------------------restart the timer */
        parms_TimerRestart(tm);

        /*--------------------Solve the linear equation */
        parms_SolverApply_b(solver, rhs, x);
        //~ int mnnz, pnnz;
        //~ parms_OperatorGetNnz(solver, &mnnz, &pnnz);
        //~ printf("mnnz %d, pnnz %d\n", mnnz, pnnz);

        //MPI_Barrier(MPI_COMM_WORLD);
        //exit(1);

        /*--------------------get total time spent on creating the pc and solving the linear
     system */
        ttol = parms_TimerGet(tm);

        printf("The time for solving  %f in proc %d\n", ttol-tpc, myid);
        printf("The total time cost  %f in proc %d\n", ttol, myid);

        /*--------------------Get the number of iterations */
        its = parms_SolverGetIts(solver);

        /*--------------------Compute the residual error  */
        parms_MatVec(A, x, y);
        for(i=0; i<llsize; i++)//nloc
        {
            y[i] = rhs[i] - y[i];
        }
        parms_VecGetNorm2_b(y, &res1, map);
        /*--------------------processor 0 outputs the result */
        if (myid == 0) {
            parms_PCGetName(pc, &name);
            parms_PCILUGetName(pc, &iluname);
            fprintf(fp, "The blocking information Bsize=%-7f, Bdensity=%-7f\n", blocksize, Bdensity);
            fprintf(fp, "The time for finding blocks=%-7f, The total blocking time=%-7f\n", tib1, tib4);
            fprintf(fp, "The preconditioner  %9s %s\n", " ", name);
            fprintf(fp, "The local preconditioner %4s %s\n", " ", iluname);
            fprintf(fp, "The memory usage %12s %-4.2f\n", "=", ratio);
            fprintf(fp, "The number of processors %4s %-4d\n", "=", npro);
            fprintf(fp, "The number of iterations %4s %-4d\n", "=", its);
            sprintf(buf, "%8.2fs", tpc);
            fprintf(fp, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol-tpc);
            fprintf(fp, "The solving time %12s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol);
            fprintf(fp, "The total time %14s %-s\n", "=", strtok(buf, " "));
            fprintf(fp, "The initial residual norm %3s %-8.2e\n", "=",norm);
            fprintf(fp, "The final residual norm %5s %-8.2e\n","=", res1);

            if(fp != stdout)
                fclose(fp);
            /* --- Write output to file ---- */
            fout = fopen("output.txt", "aw");
            fprintf(fout, "\nMatrix %d: %s \n",(mat+1), curname);
            fprintf(fout, "n = %d, nnz = %d\n", n, nnz);
            fprintf(fout, "The blocking information Bsize=%-7f, Bdensity=%-7f\n", blocksize, Bdensity);
            fprintf(fout, "The time for finding blocks=%-7f, The total blocking time=%-7f\n", tib1, tib4);
            fprintf(fout, "The preconditioner  %9s %s\n", " ", name);
            fprintf(fout, "The local preconditioner %4s %s\n", " ", iluname);
            fprintf(fout, "The memory usage %12s %-4.2f\n", "=", ratio);
            fprintf(fout, "The number of processors %4s %-4d\n", "=", npro);
            fprintf(fout, "The number of iterations %4s %-4d\n", "=", its);
            sprintf(buf, "%8.2fs", tpc);
            fprintf(fout, "The time for pc setup %7s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol-tpc);
            fprintf(fout, "The solving time %12s %-s\n", "=", strtok(buf, " "));
            sprintf(buf, "%8.2fs", ttol);
            fprintf(fout, "The total time %14s %-s\n", "=", strtok(buf, " "));
            fprintf(fout, "The initial residual norm %3s %-8.2e\n", "=",norm);
            fprintf(fout, "The final residual norm %5s %-8.2e\n","=", res1);
            fclose(fout);
        }


        //printf("after getratio\n");
        //MPI_Barrier(MPI_COMM_WORLD);
        //exit(1);
        /*
  parms_ViewerCreate(&sv, "solver.out");
  parms_TimerView(tm, sv);
*/
        /*--------------------Free memories */
        cleanCS( csmat );
        cleanVBMat( vbmat );//int cleanVBMat( vbsptr vbmat )
        free(x);
        free(y);
        free(rhs);
        free(perm);
        parms_MatFree_b(&A);
        parms_MapFree(&map);
        parms_PCFree_b(&pc);
        parms_SolverFree_b(&solver);
        parms_TimerFree(&tm);
        //printf("after cleanning");
        /*
  parms_ViewerFree(&sv);
*/
        /*----Goto next matrix ---*/
        //label1000://for block structure detection
        mat++;

    }

    /*--------------Free prm--------------*/
    free(prm);
    /*--------------------Exit MPI environment */
    MPI_Finalize();
    return 0;
}
