#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "aux.h"
#include "mmio.h"
#include <dlfcn.h>


int read_param(char *fname, char mname[MAX_MAT][MAX_LINE], fprm prm)
{
    FILE *fp, *fmat;
    char buf[MAX_LINE], *tmpstring;
    char *p, *p1, *p2;
    int i;

    if (NULL == (fp = fopen(fname, "r"))) {
        fprintf(stderr, "cannot open file %s\n", fname);
        exit(1);
    }

    for (i = 0; i < 18; i++) {
        prm->ipar[i] = -1;
    }

    /*----Get parallel Precon Type ---*/
    if(fgets(buf, MAX_LINE, fp) == NULL){
        fprintf(stderr, "Error reading global precon type");
        exit(1);
    }
    STR2U(p, p1);
    if (!strcmp(p1, "RAS")) {
        prm->pctype = PCRAS;
    }
    else if (!strcmp(p1, "SCHUR")){
        prm->pctype = PCSCHUR;
    }
    else if (!strcmp(p1, "BJ")) {
        prm->pctype = PCBJ;
    }
    free(p1);
    /* --- Get Local Precon type ---*/
    if(fgets(buf, MAX_LINE, fp) == NULL){
        fprintf(stderr, "Error reading local precon type");
        exit(1);
    }
    STR2U(p, p1);
    if (!strcmp(p1, "ILU0")) {
        prm->pcilutype = PCILU0;
    }
    else if (!strcmp(p1, "ILUK")) {
        prm->pcilutype = PCILUK;
    }
    else if (!strcmp(p1, "ILUT")) {
        prm->pcilutype = PCILUT;
    }
    else if (!strcmp(p1, "ARMS")) {
        prm->pcilutype = PCARMS;
    }
    else if (!strcmp(p1, "VBARMS")) {
        prm->pcilutype = PCVBARMS;
    }
    else if (!strcmp(p1, "VBARMSOLD")) {
        prm->pcilutype = PCVBARMSOLD;
    }
    else if (!strcmp(p1, "VBILUT")) {
        prm->pcilutype = PCVBILUT;
    }
    else if (!strcmp(p1, "VBILUK")) {
        prm->pcilutype = PCVBILUK;
    }
    free(p1);
    /*--- Get tolerance for local solve---*/
    if(fscanf(fp,"%lf",&prm->pgfpar[0]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Get tolerance for global solve ---*/
    if(fscanf(fp,"%lf",&prm->pgfpar[1]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Get number of levels ---*/
    if(fscanf(fp,"%d",&prm->ipar[0]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Use symmetric or nonsymmetric perm ---*/
    if(fscanf(fp,"%d",&prm->ipar[1]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- B block size ---*/
    if(fscanf(fp,"%d",&prm->ipar[2]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Tolerance for independent sets ---*/
    if(fscanf(fp,"%lf",&prm->tolind) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Outer Krylov dimension ---*/
    if(fscanf(fp,"%d",&prm->ipar[6]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Maximum outer iterations ---*/
    if(fscanf(fp,"%d",&prm->ipar[7]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Inner Krylov dimension ---*/
    if(fscanf(fp,"%d",&prm->ipar[4]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- maximum number of inner iterations ---*/
    if(fscanf(fp,"%d",&prm->ipar[5]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- nonsymmetric perm. for interlevel blocks ---*/
    if(fscanf(fp,"%d",&prm->ipar[10]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Column Permutations - interlevel blocks ---*/
    if(fscanf(fp,"%d",&prm->ipar[11]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Row scaling - interlevel blocks ---*/
    if(fscanf(fp,"%d",&prm->ipar[12]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Column scaling - interlevel blocks ---*/
    if(fscanf(fp,"%d",&prm->ipar[13]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Nonsymm. Perm. for last level SC. ---*/
    if(fscanf(fp,"%d",&prm->ipar[14]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Col. Perm. for last level SC. ---*/
    if(fscanf(fp,"%d",&prm->ipar[15]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Row scaling - last level SC. ---*/
    if(fscanf(fp,"%d",&prm->ipar[16]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- Column scaling - last level SC. ---*/
    if(fscanf(fp,"%d",&prm->ipar[17]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- lfil0 for ilut, iluk, and arms for lfil[0-3] ---*/
    if(fscanf(fp,"%d",&prm->lfil[0]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- lfil for shur ---*/
    if(fscanf(fp,"%d",&prm->lfil[4]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- lfil for ILUT L and U ---*/
    if(fscanf(fp,"%d",&prm->lfil[5]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- droptol0(droptol[0=3], L, U, L^{-1}F, EU^{-1}---*/
    if(fscanf(fp,"%lf",&prm->droptol[0]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- droptol4(for schur complements at each level) ---*/
    if(fscanf(fp,"%lf",&prm->droptol[4]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- droptol5(for ILUT in last level Schur Complement) ---*/
    if(fscanf(fp,"%lf",&prm->droptol[5]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- eps(For block version local solver(vbarms and vbilut)) ---*/
    if(fscanf(fp,"%lf",&prm->eps) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- for the cosine graph compression algorithm ---*/
    if(fscanf(fp,"%d",&prm->cosine) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- for the cosine graph compression algorithm ---*/
    if(fscanf(fp,"%d",&prm->lfil[7]) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    while(fgetc(fp) != '\n');
    /*--- for a constant block size ---*/
    if(fscanf(fp,"%d",&prm->constant_block_size) != 1){
        fprintf(stderr, "Error reading value");
        exit(1);
    }
    prm->droptol[1] = prm->droptol[2] = prm->droptol[3] =
            prm->droptol[0];
    prm->droptol[6] = prm->droptol[5];

    prm->lfil[1] = prm->lfil[2] = prm->lfil[3] =
            prm->lfil[0];
    prm->lfil[6] = prm->lfil[5];

    fclose(fp);
    
    /*------- Done reading inputs, now read list of Matrices -------------*/
#if defined(DBL_CMPLX)
    if (NULL == (fmat = fopen("./matfileCmplx", "r"))) {
        fprintf(stderr, "cannot open file: matfile\n");
        exit(1);
    }
#else  
    if (NULL == (fmat = fopen("./matfileReal", "r"))) {
        fprintf(stderr, "cannot open file: matfile\n");
        exit(1);
    }
#endif  
    i = 0;
    do {
        if(fgets(buf,MAX_LINE,fmat) == NULL){
            fprintf(stderr, "Error reading matrix list");
            exit(1);
        }
        tmpstring = (char *) strtok(buf, "\n");
        if (tmpstring == NULL)
            continue;
        strcpy(mname[i], tmpstring);
        i++;
    } while (strncmp(buf, "##", 2));

    fclose(fmat);

    return 0;
}


void set_pc_params_b(parms_PC pc, fprm prm)//block version
{
    parms_PCSetType_b(pc,          prm->pctype);
    parms_PCSetILUType(pc,         prm->pcilutype);
    parms_PCSetNlevels(pc,         prm->ipar[0]);
    parms_PCSetPermType(pc,        prm->ipar[1]);
    parms_PCSetBsize(pc,           prm->ipar[2]);
    parms_PCSetInnerEps(pc,        prm->pgfpar[0]);
    parms_PCSetInnerKSize(pc,      prm->ipar[4]);
    parms_PCSetInnerMaxits(pc,     prm->ipar[5]);
    parms_PCSetFill(pc,            prm->lfil);
    parms_PCSetTol(pc,             prm->droptol);
    parms_PCSetTolInd(pc,          prm->tolind);
    parms_PCSetPermScalOptions(pc, &prm->ipar[10], 1);
    parms_PCSetPermScalOptions(pc, &prm->ipar[14], 0);
}


void set_solver_params(parms_Solver solver, fprm prm)
{
    char buf[BUFFLEN];

    //printf("ipar[7] = %d\n", prm->ipar[7]);
    sprintf(buf, "%d", prm->ipar[7]);
    parms_SolverSetParam(solver, MAXITS, buf);
    //printf("ipar[6] = %d\n", prm->ipar[6]);
    sprintf(buf, "%d", prm->ipar[6]);
    parms_SolverSetParam(solver, KSIZE,  buf);
    sprintf(buf, "%g", prm->pgfpar[1]);
    parms_SolverSetParam(solver, DTOL,  buf);
}

void fread_param_(char *fname, fprm *prm, char *matrix, int *matlen, int len)
{
    char mname[MAX_MAT][MAX_LINE], *buff, *buff2;

    buff2 = malloc((len+1)*sizeof(*buff2));
    strncpy(buff2, fname, len);
    buff2[len] = '\0';
    *prm = malloc(sizeof(**prm));
    read_param(buff2, mname, *prm);
    buff = &mname[0][0];
    strncpy(matrix, buff, strlen(buff));
    matrix[strlen(buff)] = '\0';
    *matlen = (int)strlen(matrix);
    free(buff2);

}


void fset_solver_params_(parms_Solver *solver, fprm *prm)
{
    set_solver_params(*solver, *prm);
}

void fprm_free_(fprm *prm)
{
    free(*prm);
}

//-------------------------------------dividing line, the subroutines below are added by George Liao-----------------------------------------------------------------------------------------------------
void coo2csptr(int n, int nnz, FLOAT *a, int *ir, int *jc, csptr mat){

    /* ----------------------------------------------------------------------- */
    /*   Coordinate     to   Compressed Sparse Row C-style */
    /* -----------------------------------------------------------------------  */
    /*  converts a matrix that is stored in coordinate format */
    /*   a, ir, jc into a row general sparse csptr format. */

    /*  on entry: */
    /* ---------  */
    /*  n	= dimension of the matrix  */
    /*  nnz	= number of nonzero elements in matrix */
    /*  a, */
    /*  ir,  */
    /*  jc    = matrix in coordinate format. a(k), ir(k), jc(k) store the nnz */
    /*          nonzero elements of the matrix with a(k) = actual real value of */
    /*  	  the elements, ir(k) = its row number and jc(k) = its column  */
    /* 	  number. The order of the elements is arbitrary. */

    /* on return: */
    /* -----------  */
    /*  ir 	is destroyed */
    /*  */
    /*  csptr mat = matrix in general sparse C-style matrix format*/
    int *counter, len, i, *pja;
    FLOAT *paa;


    counter = malloc(n*sizeof(*counter));
    memset(counter, 0, n*sizeof(int));

    setupCS(mat,n,1);//int setupCS(csptr amat, int len, int job)

    memset(mat->nnzrow, 0, n*sizeof(int));
    
    for (i = 0; i < nnz; i++)
        mat->nnzrow[ir[i]-1]++;
    
    for(i = 0; i < n; i++){
        len = mat->nnzrow[i];
        if (len) {
            pja = malloc(len*sizeof(*pja));
            paa = malloc(len*sizeof(*paa));
            mat->pj[i] = pja;
            mat->pa[i] = paa;
        }
    }

    for (i = 0; i < nnz; i++){
        mat->pj[ir[i]-1][counter[ir[i]-1]] = jc[i]-1;
        mat->pa[ir[i]-1][counter[ir[i]-1]] = a[i];
        counter[ir[i]-1]++;
    }

    free(counter);counter = NULL;

}



void output_intvector(char *filename,int *v,int i0, int i1){
    //output a integer vector
    FILE *fp;
    int jj;
    //    printf(" in side value is.\n" );//%f %p %s %c

    fp = fopen(filename,"w");
    for(jj=i0; jj<i1;jj++)
        fprintf(fp,"%d\n",v[jj]);
    fclose(fp);

}

void output_csrmatrix(char *filename,int *ia, int *ja, double *a, int n){

    FILE *fp;
    int jj, nnz;

    fp = fopen(filename, "w");

    fprintf(fp,"%d\n", n);

    for(jj = 0; jj < n+1; jj++)
        fprintf(fp, "%d\n", ia[jj]);

    nnz = ia[n] - 1;

    printf("nnz value is %d.\n", nnz);//%f %p %s %c

    for(jj = 0; jj < nnz; jj++)
        fprintf(fp, "%d\n", ja[jj]);

    for(jj = 0; jj < nnz; jj++)
        fprintf(fp, "%20.16e\n", a[jj]);

    fclose(fp);

}

void output_dblvector(char *filename, FLOAT *v, int i0, int i1){//output a double vector 
    int jj;
    FILE *fp=fopen(filename,"w");

#if defined(DBL_CMPLX)
    for(jj=i0; jj<i1;jj++)
        fprintf(fp,"%20.16e  %20.16e\n",creal(v[jj]),cimag(v[jj]));
#else
    for(jj=i0; jj<i1;jj++)
        fprintf(fp, "%20.16e\n", v[jj]);//%20.16e//%f\n

#endif
    fclose(fp);

}

int outputcsmat ( csptr mat, char *filename, int onebase){//new
    /*----------------------------------------------------------------------
| Output the pattern of CSR format matrix, which can be loaded by matlab
  input: onebase means the coordinate in ouput coofile start with 0 or 1.
----------------------------------------------------------------------*/
    FILE *fmatlab = fopen( filename, "w" );
    int n = mat->n, i, j, nnzrow, pid;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if ((onebase!=0) &&(onebase!=1)){
        printf("input parameter error\n");
        return 1;
    }
    printf("\n outputcsmat n=%d\n",n);//getchar();
    if( !fmatlab ) {
        printf("error outputcsmat"); getchar();
        return -1;
    }
    fprintf( stdout, "%d %d %d pid=%d\n", n, n, nnzCS(mat),pid );
    fprintf( fmatlab, "%d %d %d\n", n, n, nnzCS(mat) );//getchar();
    /*old:
  fprintf( stdout, "%d %d %d\n", n, n, cs_nnz(mat) );
  fprintf( fmatlab, "%d %d %d\n", n, n, cs_nnz(mat) );getchar();
*/


    for( i = 0; i < n; i++ ) {
        //fprintf(stdout,"i=%d"); if((i+1)% 5==0)printf("\n");
        nnzrow = mat->nnzrow[i];
        for( j = 0; j < nnzrow; j++ )
#if defined(DBL_CMPLX)
            fprintf( fmatlab, "%d %d %20.16e %20.16e\n", i+onebase, mat->pj[i][j]+onebase, creal(mat->pa[i][j]),cimag(mat->pa[i][j]) );
#else
            //if ((mat->pj[i][j]+1)>0)
            fprintf( fmatlab, "%d %d %20.16e\n", i+onebase, mat->pj[i][j]+onebase,mat->pa[i][j] );
#endif

    }

    fclose( fmatlab );
    return 0;
}



int colunms2csptr(int n, int *ia, int *ja, FLOAT *a, csptr mat)
{
    /*
 * This subroutine converts the 3 colunms CSR format to C-style CSR format
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 *
 * Three arrars: ia, ja, a.
 *
 * on return:
 *-----------
 *
 *  mat = Sparse Row format Matrix, csptr
 *
 *  ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int i, j, firstpos;
    int *nnzrow, *pj;

    setupCS(mat,n,1);//int setupCS(csptr amat, int len, int job)

    for( i = 0; i < n; i++ )
    {
        nnzrow = mat->nnzrow;
        nnzrow[i] = B_DIM(ia,i);//get the lenth of each row
        firstpos = ia[i]-1;//get the address of the first entry in each row
        mat->pj[i] = (int *) Malloc(nnzrow[i]*sizeof(int), "colunms2csptr" );//ja=mat->ja[i]=vbmat->ja[i];
        pj = mat->pj[i];
        for( j = 0; j < nnzrow[i]; j++ )
        {
            pj[j] = ja[firstpos+j]-1;//ja[firstpos+j]--;
        }
        mat->pa[i] = (FLOAT *) Malloc(nnzrow[i]*sizeof(FLOAT), "colunms2csptr" );//ja=mat->ja[i]=vbmat->ja[i];
        memcpy(mat->pa[i],&a[firstpos],nnzrow[i]*sizeof(FLOAT));
    }
    return 0;
}


int csptr2colunms(csptr mat, int *n, int *ia, int *ja)//
{
    /*
 * This subroutine extracts the index vectors from the C-style CSR matrix.
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 *
 *  mat = Sparse Row format Matrix, csptr
 *
 * on return:
 *-----------
 *
 *  Two arrays and one interger: ia, ja, n.
 *
 *  ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int i, j, firstpos;
    int *nnzrow, *pj;


    ia[0] = 1;
    *n = mat->n;
    for( i = 0; i < mat->n; i++ )
    {
        nnzrow = mat->nnzrow;
        ia[i+1] = ia[i] + nnzrow[i];
        firstpos = ia[i]-1;//get the address of the first entry in each row
        pj = mat->pj[i];
        for( j = 0; j < nnzrow[i]; j++ )
        {
            ja[firstpos+j] = pj[j] + 1;//pj[j] = ja[firstpos+j]-1;//ja[firstpos+j]--;
        }
    }
    return 0;
}

int outputvbmat1( vbsptr vbmat, char *filename, int onebase){
    /*----------------------------------------------------------------------
| Output the pattern of vbsptr, which can be loaded by matlab for rectangular block matrix only
----------------------------------------------------------------------*/
    FILE *fmatlab = fopen( filename, "w" );
    int i, j, nzcount, col, dimR, dimC ,k, pid;
    int bn = vbmat->n, bnc = vbmat->nc, *ja, *bsz = vbmat->bsz, *bszc=vbmat->bszc, n=bsz[bn], nc=bszc[bnc];
    BData *ba;

    MPI_Comm_rank(MPI_COMM_WORLD, &pid);

    if ((onebase!=0) &&(onebase!=1)){
        printf("input parameter error\n");
        return 1;
    }
    printf("\n outputvbmat1 n=%d\n",n);//getchar();
    if( !fmatlab ) {
        printf("error outputvbmat1"); getchar();
        return -1;
    }
    //fprintf( stdout, "%d %d %d\n", n, nc, memVBMat1(vbmat) );
    fprintf( stdout, "%d %d %d pid=%d\n", n, nc, memVBMat1(vbmat),pid );
    fprintf( fmatlab, "%d %d %d\n", n, nc, memVBMat1(vbmat) );getchar();
    /*old:
  fprintf( stdout, "%d %d %d\n", n, n, cs_nnz(mat) );
  fprintf( fmatlab, "%d %d %d\n", n, n, cs_nnz(mat) );getchar();
*/
    printf("total nonzeroblocks = %d, vbmat->lsize = %d\n",nnzVBMat(vbmat),vbmat->lsize); getchar();

    for( i = 0; i < bn; i++ ) {
        dimR = B_DIM(bsz,i);//get the row dimension of the block
        nzcount = vbmat->nzcount[i];//
        ja = vbmat->ja[i];
        ba = vbmat->ba[i];
        //fprintf(stdout,"i=%d"); if((i+1)% 5==0)printf("\n")
        for( j = 0; j < nzcount; j++ ){
            //printf("\n i=%d,j=%d, nzcount=%d,col=%d",i,j,nzcount,col); getchar();
            col = ja[j];
            //printf("\n col=%d",col); getchar();
            dimC = B_DIM(bszc,col-vbmat->lsize);//obtain the column dimension of the block
            for( k = 0; k < dimC*dimR; k++ )
#if defined(DBL_CMPLX)
                fprintf( fmatlab, "%d %d %20.16e %20.16e\n", bsz[i]+k%dimR+onebase, bszc[col-vbmat->lsize]+k/dimR+onebase, creal(ba[j][k]),cimag(ba[j][k]) );
#else
                fprintf( fmatlab, "%d %d %20.16e\n", bsz[i]+k%dimR+onebase, bszc[col-vbmat->lsize]+k/dimR+onebase,ba[j][k] );
#endif
        }
    }

    fclose( fmatlab );
    return 0;
}

int outputvbmat( vbsptr vbmat, char *filename, int onebase){
    /*----------------------------------------------------------------------
| Output the pattern of vbsptr matrix, which can be loaded by matlab
----------------------------------------------------------------------*/

    if(vbmat->bszc==NULL){

        FILE *fmatlab = fopen( filename, "w" );
        int i, j, nzcount, col, dimR, dimC ,k, pid;
        int bn = vbmat->n, *ja, *bsz = vbmat->bsz,n=bsz[bn];
        BData *ba;

        MPI_Comm_rank(MPI_COMM_WORLD, &pid);

        if ((onebase!=0) &&(onebase!=1)){
            printf("input parameter error\n");
            return 1;
        }
        printf("\n outputvbmat n=%d\n",n);//getchar();
        if( !fmatlab ) {
            printf("error outputvbmat"); getchar();
            return -1;
        }
        fprintf( stdout, "%d %d %d pid=%d\n", n, n, memVBMat(vbmat),pid );
        fprintf( fmatlab, "%d %d %d\n", n, n, memVBMat(vbmat) );//getchar();
        /*old:
    fprintf( stdout, "%d %d %d\n", n, n, cs_nnz(mat) );
    fprintf( fmatlab, "%d %d %d\n", n, n, cs_nnz(mat) );getchar();
    */
        //printf("%d\n",nnzVBMat(vbmat)); getchar();
        printf("total nonzeroblocks = %d, vbmat->lsize = %d\n",nnzVBMat(vbmat),vbmat->lsize); //getchar();

        for( i = 0; i < bn; i++ ) {
            dimR = B_DIM(bsz,i);//get the row dimension of the block
            nzcount = vbmat->nzcount[i];//
            ja = vbmat->ja[i];
            ba = vbmat->ba[i];
            //fprintf(stdout,"i=%d"); if((i+1)% 5==0)printf("\n")
            for( j = 0; j < nzcount; j++ ){
                col = ja[j];
                // printf("\n i=%d,j=%d, nzcount=%d,col=%d",i,j,nzcount,col); getchar();
                //printf("\n col=%d",col); getchar();
                dimC = B_DIM(bsz,col-vbmat->lsize);//obtain the column dimension of the block
                for( k = 0; k < dimC*dimR; k++ ){
#if defined(DBL_CMPLX)
                    fprintf( fmatlab, "%d %d %20.16e %20.16e\n", bsz[i]+k%dimR+onebase, bsz[col-vbmat->lsize]+k/dimR+onebase, creal(ba[j][k]),cimag(ba[j][k]) );
                    //if(fabs(ba[j][k]))kk++;
#else
                    fprintf( fmatlab, "%d %d %20.16e\n", bsz[i]+k%dimR+onebase, bsz[col-vbmat->lsize]+k/dimR+onebase,ba[j][k] );
#endif
                }
                //printf("kk=%d",kk);getchar();
                //if(kk==0)printf("zero block encountered");
            }
        }

        fclose( fmatlab );
        return 0;
    }
    else return outputvbmat1(vbmat, filename, onebase);
}


int vbsptr2colunms(vbsptr mat, int *bnnz, int *ia, int *ja, BData *a)//
{
    /*
 * This subroutine converts C-style VBCSR matrix to 3 columns format matrix.
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 *
 *  mat = Sparse Row format Matrix, csptr
 *
 * on return:
 *-----------
 *
 *  Three arrays and one interger: ia, ja, a, bnnz.
 *
 *  ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int i, j, firstpos;
    int *nzcount, *pj;
    BData *pa;

    ia[0] = 1;

    for( i = 0; i < mat->n; i++ )
    {

        nzcount = mat->nzcount;
        ia[i+1] = ia[i] + nzcount[i];
        firstpos = ia[i]-1;//get the address of the first entry in each row
        pj = mat->ja[i];
        pa = mat->ba[i];
        for( j = 0; j < nzcount[i]; j++ )
        {
            //---------------------------------------------------------------------------------------------------------
            ja[firstpos+j] = pj[j] + 1;
            a[firstpos+j] = pa[j];

        }
    }
    *bnnz = ia[mat->n]-1;
    return 0;
}




int outputcsmatpa ( csptr mat, char *filename, int onebase){//for parallel csr matrix output

    int length, myid;
    char *newfilename;//*filename = "A->aux_data",
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    length = strlen(filename) + 6;
    length += myid / 10 + 1;
    PARMS_NEWARRAY(newfilename, length);

    strcpy(newfilename, filename);
    sprintf(newfilename, "%s_%d.coo", newfilename, myid);
    outputcsmat(mat,newfilename,onebase);//getchar();
    PARMS_FREE(newfilename);
    if (filename == NULL) {
        return -1;
    }
    else {
        return 0;
    }
}

int output_intvectorpa(char *filename,int *v,int i0, int i1){//for parallel integer vector output

    int length, myid;
    char *newfilename;//*filename = "A->aux_data",
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    length = strlen(filename) + 6;
    length += myid / 10 + 1;
    PARMS_NEWARRAY(newfilename, length);

    strcpy(newfilename, filename);
    sprintf(newfilename, "%s_%d.coo", newfilename, myid);
    output_intvector(newfilename,v,i0,i1);//getchar();
    PARMS_FREE(newfilename);
    if (filename == NULL) {
        return -1;
    }
    else {
        return 0;
    }
}

int output_dblvectorpa(char *filename,FLOAT *v,int i0, int i1){//for parallel integer vector output

    int length, myid;
    char *newfilename;//*filename = "A->aux_data",
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    length = strlen(filename) + 6;
    length += myid / 10 + 1;
    PARMS_NEWARRAY(newfilename, length);

    strcpy(newfilename, filename);
    sprintf(newfilename, "%s_%d.coo", newfilename, myid);
    output_dblvector(newfilename,v,i0,i1);//getchar();
    PARMS_FREE(newfilename);
    if (filename == NULL) {
        return -1;
    }
    else {
        return 0;
    }
}

int outputvbmatpa ( vbsptr mat, char *filename, int onebase){//for parallel vbcsr matrix output

    int length, myid;
    char *newfilename;//*filename = "A->aux_data",
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    length = strlen(filename) + 6;
    length += myid / 10 + 1;
    PARMS_NEWARRAY(newfilename, length);

    strcpy(newfilename, filename);
    //printf("A->aux_data=%d",A->aux_data);getchar();
    sprintf(newfilename, "%s_%d.coo", newfilename, myid);
    outputvbmat(mat,newfilename,onebase);//getchar();
    PARMS_FREE(newfilename);
    if (filename == NULL) {
        return -1;
    }
    else {
        return 0;
    }
}


//---------------------------------------------------binary file reading---------------------------------------------------------------

int bincols2csptr(int n, int *nnzptr, int *ja, FLOAT *a, csptr mat)
{
    /*
 * This subroutine converts the 3 colunms CSR format to C-style CSR format
 *
 *----------------------------------------------------------------------
 * on entry:
 *----------
 *
 * Three arrays: nnzptr, ja, a.
 * nnzptr is a array to store the number of non-zero entries on each row.
 *
 * on return:
 *-----------
 *
 *  mat = Sparse Row format Matrix, csptr
 *
 *  ierr  = integer, error code.
 *              0  -- normal termination
 *             -1  -- error occur
 *
 *---------------------------------------------------------------------*/
    int i;
    int firstpos = 0;
    setupCS(mat, n ,1);//int setupCS(csptr amat, int len, int job)
    int *nnzrow = mat->nnzrow;

    for (i = 0; i < n; ++i)
    {
        nnzrow[i] = nnzptr[i];//get the lenth of each row
        mat->pj[i] = (int *) Malloc(nnzrow[i]*sizeof(int), "colunms2csptr" );//ja=mat->ja[i]=vbmat->ja[i];
        memcpy(mat->pj[i], &ja[firstpos], nnzrow[i]*sizeof(int));
        mat->pa[i] = (FLOAT*) Malloc(nnzrow[i]*sizeof(FLOAT), "colunms2csptr" );//ja=mat->ja[i]=vbmat->ja[i];
        memcpy(mat->pa[i],&a[firstpos],nnzrow[i]*sizeof(FLOAT));
        firstpos += nnzrow[i];
    }
    return 0;
}


double bswap_double(double a)//
{
    double b;
    unsigned char *src = (unsigned char *)&a,
            *dst = (unsigned char *)&b;

    int i;
    for (i = 0; i < 8; ++i) {
        dst[i] = src[7-i];
    }
    return b;

}

double byteswap_double(double v) // This doesn't work for some values
{
    union { // This trick is first used
        int64_t i;
        double  d;
    } conv;
    conv.d = v;
    conv.i = bswap_64(conv.i);
    return conv.d;
}


int fillinuppertrg(int n, int nnz, FLOAT *a, int *ja, int *ia, FLOAT *b, int *jb, int *ib, SYMMTYPE sytype){
    //fill in the upper triangle part when matrix is SYMM or HERMITION or SKEW_SYMM
    int i, j, col;
    int *nzcount = malloc(n*sizeof(*nzcount));

    for(i = 0; i < n; i++)
        nzcount[i] = ib[i+1] - ib[i];

    for(i = 0; i < n; i++) {
        for(j = ib[i]-1; j < ib[i+1]-1; j++) {
            col = jb[j] - 1;
            if(col != i) nzcount[col]++;
        }
    }

    for(i = 0; i < n; i++)
        ia[i+1] = ia[i] + nzcount[i];
    //symmetric

    memset(nzcount, 0, n * sizeof(int));
    switch (sytype) {
    case SYMM:
        for(i = 0; i < n; i++) {
            for(j = ib[i]-1; j < ib[i+1]-1; j++){
                col = jb[j] - 1;
                ja[ia[i]-1+nzcount[i]] = jb[j];
                a[ia[i]-1+nzcount[i]] = b[j];
                nzcount[i]++;
                if(col != i) {
                    ja[ia[col]-1+nzcount[col]] = i+1;
                    a[ia[col]-1+nzcount[col]] = b[j];
                    nzcount[col]++;
                }
            }
        }
        break;
    case SKEW_SYMM:
        for(i = 0; i < n; i++) {
            for(j = ib[i]-1; j < ib[i+1]-1; j++){
                col = jb[j] - 1;
                ja[ia[i]-1+nzcount[i]] = jb[j];
                a[ia[i]-1+nzcount[i]] = b[j];
                nzcount[i]++;
                if(col != i) {
                    ja[ia[col]-1+nzcount[col]] = i+1;
                    a[ia[col]-1+nzcount[col]] = - b[j];
                    nzcount[col]++;
                }
            }
        }
        break;
    case HERMITION:
#if defined(DBL_CMPLX)
        for(i = 0; i < n; i++) {
            for(j = ib[i]-1; j < ib[i+1]-1; j++){
                col = jb[j] - 1;
                ja[ia[i]-1+nzcount[i]] = jb[j];
                a[ia[i]-1+nzcount[i]] = b[j];
                nzcount[i]++;
                if(col != i) {
                    ja[ia[col]-1+nzcount[col]] = i+1;
                    a[ia[col]-1+nzcount[col]] = conj(b[j]);
                    nzcount[col]++;
                }
            }
        }
#endif
        break;
    default:
        break;
    }
    free(nzcount);
    return 0;
}


int fullmatize(int n, int nnz, FLOAT *a, int *ja, int *ia, FLOAT *b, int *jb, int*ib, char *type){


    if(type[0] == 'R' || type[0] == 'r'){
        switch (type[1]) {
        case 'S':
        case 's':
            //symmetric
            fillinuppertrg(n, nnz, a, ja, ia, b, jb, ib, SYMM);
            break;
            //skew symmetric
        case 'Z':
        case 'z':
            fillinuppertrg(n, nnz, a, ja, ia, b, jb, ib, SKEW_SYMM);
            break;
        case 'U':
        case 'u':
            //unsymmetric
            memcpy(ia, ib, (n+1)*sizeof(*ia));
            memcpy(ja, jb, nnz*sizeof(*ja));
            memcpy(a,  b, nnz*sizeof(*a));
            break;
        case 'R':
        case 'r':
            //rectangular matrix
            fprintf(stderr, "Rectangular matrix, can not not process\n");
            return 1;
        default:
            fprintf(stderr, "Incorrect matrix type\n");
            return 1;
            //   break;
        }
    }

    else if(type[0] == 'C' || type[0] == 'c'){
        switch (type[1]) {
        case 'S':
        case 's':
            //symmetric
            printf("in symm case \n");//%f %p %s %c

            fillinuppertrg(n, nnz, a, ja, ia, b, jb, ib, SYMM);
            break;
        case 'H':
        case 'h':
            //Hermition
            printf("in Hermition case \n");//%f %p %s %c
            fillinuppertrg(n, nnz, a, ja, ia, b, jb, ib, HERMITION);
            break;
        case 'U':
        case 'u':
            //unsymmetric
            memcpy(ia, ib, (n+1)*sizeof(*ia));
            memcpy(ja, jb, nnz*sizeof(*ja));
            memcpy(a,  b, nnz*sizeof(*a));
            break;
        case 'R':
        case 'r':
            //rectangular matrix
            fprintf(stderr, "Rectangular matrix, can not process\n");
            return 1;
        default:
            fprintf(stderr, "Incorrect matrix type\n");
            return 1;
            //   break;
        }
    }
    else
        fprintf(stderr, "Neither real or complex matrix, can not process\n");

    //    free(nzcount);
    return 0;
}

int compare(const void * a, const void * b)//for qsort
{
    return ( *(int*)a - *(int*)b );
}


int local_read_bin_data_from_indices(FILE *binfile, long int M, long int N, int nz, int nnzptr[], int ja[],
                                     FLOAT val[], int *indices, int nindices, int *gia)
{

    long long offset, ja_startindex, gpindex, numread;

    int i1, rowlength;
    ja_startindex = 0;

    memset(nnzptr, 0, sizeof(int)*M);

    for (i1 = 0; i1 < nindices; i1++) {
        gpindex = indices[i1];
        rowlength = nnzptr[gpindex] = gia[gpindex+1] - gia[gpindex];

        offset = (4 + M + gia[gpindex])*sizeof(int);//move pointer to load ja's one row colunm index
        fseek(binfile, offset, 0);
        numread = fread(&ja[ja_startindex], sizeof(int), rowlength, binfile);
        if (numread != rowlength){
            fprintf(stderr," error in fread");
            MPI_Finalize();
            exit(1);
        }

        offset = (4 + M + gia[M])*sizeof(int) + gia[gpindex]*sizeof(FLOAT);//move pointer to load ja's one row colunm index
        fseek(binfile, offset, 0);

        numread = fread(&val[ja_startindex], sizeof(FLOAT), rowlength, binfile);
        if (numread != rowlength){
            fprintf(stderr," error in fread");
            MPI_Finalize();
            exit(1);
        }
        ja_startindex += rowlength;
    }

    for (i1 = 0; i1 < ja_startindex; i1++){
        ja[i1] = bswap_32(ja[i1]);//swap between big endian and small endian
        val[i1] = byteswap_double(val[i1]);
    }

    return 0;

}

int local_read_bin_data_b(FILE *binfile, long int M, long int N, int nz, int nnzptr[], int ja[],
                          FLOAT val[], int *idom, int *dom, int *perm, int *nB, int nBlock, int *gia)
{
    int *bsz = (int *) calloc(nBlock+1, sizeof(int));

    int i1, gpindex, gbindex;
    int i2 = 0;
    int myid;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    for (i1 = 0; i1 < nBlock; ++i1)
        bsz[i1+1] = bsz[i1] + nB[i1];
    int istart = idom[myid];
    int iend = idom[myid+1];

    int *iperm = malloc(M*sizeof(*iperm));
    for (i1 = 0; i1 < M; i1++){
        iperm[perm[i1]] = i1;//generate inverse permutation array
    }

    int *indices = malloc(M*sizeof(int));
    int counter = 0;//for nnzrow index

    for (i1 = istart; i1 < iend; i1++) {
        gbindex = dom[i1-1] - 1;// global block-wise index

        for (i2 = bsz[gbindex]; i2 < bsz[gbindex+1]; i2++) {//i1 is the glocal point-wise index on block matrix
            gpindex = iperm[i2];//i2 map to gpindex
            indices[counter++] = gpindex;
        }
    }
    qsort(indices, counter, sizeof(int), compare);

    int ierr = local_read_bin_data_from_indices(binfile, M, N, nz, nnzptr, ja, val, indices, counter, gia);

    free(bsz);
    free(iperm);
    free(indices);

    return ierr;
}

int local_read_bin_data(FILE *binfile, long int M, long int N, int nz, int ia[], int ja[],
                        FLOAT val[], int *idom, int *dom, int *gia)
{
    int i1, gpindex;
    int myid;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    int istart = idom[myid];
    int iend = idom[myid+1];

    int *indices = malloc(M*sizeof(int));
    int counter = 0;//for nnzrow index

    for (i1 = istart; i1 < iend; i1++) {
        gpindex = dom[i1-1] - 1;//i2 map to gpindex
        indices[counter++] = gpindex;
    }
    qsort(indices, counter, sizeof(int), compare);

    int *nnzptr = calloc(M, sizeof(int));

    int ierr = local_read_bin_data_from_indices(binfile, M, N, nz, nnzptr, ja, val, indices, counter, gia);

    // Create ia and ja which start counting at 1. Supermooi!
    ia[0] = 1;
    for (i1 = 0; i1 < M; i1++)
        ia[i1+1] = ia[i1] + nnzptr[i1];

    for (i1 = 0; i1 < ia[M]; i1++)
        ja[i1] += 1;

    free(nnzptr);
    free(indices);

    return ierr;
}


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
