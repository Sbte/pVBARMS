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
#include<byteswap.h>//biswap_32


//double bswap_double(double a)
//{
//    //    uint64_t b;
//    double b;
//    unsigned char *src = (unsigned char *)&a,
//            *dst = (unsigned char *)&b;
//    //    if (is_littleEndian())
//    //    {
//    int i;
//    for (i = 0; i < 8; ++i) {
//        dst[i] = src[7-i];
//    }
//    //        dst[0] = src[7];
//    //        dst[1] = src[6];
//    //        dst[2] = src[5];
//    //        dst[3] = src[4];
//    //        dst[4] = src[3];
//    //        dst[5] = src[2];
//    //        dst[6] = src[1];
//    //        dst[7] = src[0];
//    //    }
//    //    else
//    //        b = *(uint64_t *)&a;
//    return b;

//}


int main(int argc, char *argv[])
{
    //    int counter;
    FILE *ptr_myfile;
    //    struct rec my_record;
    //    char buffer[20];
    long int header;
    int i;//for loop

    MPI_Init(&argc, &argv);//initilization


    ptr_myfile = fopen("mat000.bin","rb");//r and rb is for bin file
    if (!ptr_myfile)
    {
        fprintf(stderr," error in fopen");
        return 1;
    }

    printf("open file successfully!\n");

    printf("sizeof int is %ld \n", sizeof(int) );
    printf("sizeof long int is %ld \n", sizeof(long int) );
    printf("sizeof float is %ld \n", sizeof(float) );
    //    printf("sizeof uint64_t is %ld \n", sizeof(uint64_t) );
    printf("sizeof 4*int is %ld \n", sizeof(int)*4 );
    //printf("sizeof uint32_t is %ld \n", sizeof(uint32_t) );


    printf("sizeof double is %ld \n", sizeof(double) );
    printf("sizeof long double is %ld \n", sizeof(long double) );
    //				ia = malloc(nnz*sizeof(*ia));





    //    nnz = fread(fd,m,'int32');  %nonzeros per row
    //    sum_nz = sum(nnz);
    //    long int m =
    long int numread;

    numread = fread(&header, sizeof(int), 1, ptr_myfile);
    header = bswap_32(header);

    printf("header is %ld \n", header );
//    exit(1);
    if (header == 1211216){
        long int m ;    // = header(1); //row;
        long int n ;  //   = header(2); //column;
        long int nz;//     = header(3);%nnz


        numread = fread(&m, sizeof(int), 1, ptr_myfile);
        m = bswap_32(m);//row
        printf("m is %ld \n", m );

        numread = fread(&n, sizeof(int), 1, ptr_myfile);
        n = bswap_32(n);//column
        printf("n is %ld \n", n );

        numread = fread(&nz, sizeof(int), 1, ptr_myfile);
        nz = bswap_32(nz);// nnz
        printf("nz is %ld \n", nz);

        int *nnzptr = malloc(m*sizeof(*nnzptr)); //nnz of each row
        //    int *nnz = (int*) malloc (m*sizeof(int));

        numread = fread(nnzptr, sizeof(int), m, ptr_myfile);
        if (numread != m){
            fprintf(stderr," error in fread");
            return 1;
        }

        long int sum_nz = 0;
        for (i = 0; i < m; ++i){
            nnzptr[i] = bswap_32(nnzptr[i]);
            sum_nz += nnzptr[i];
            //        printf("nnz[i] = %d\n", nnz[i]);//%f %p %s %c
        }
        //    printf("sum_nz is %ld\n", sum_nz);//%f %p %s %c

        if(sum_nz != nz){
            fprintf(stderr," No-Nonzeros sum-rowlengths do not match %ld %ld", nz, sum_nz);
            return 1;
        }

        int *ja = malloc(nz*sizeof(*ja)); //column indeces of all non-zero entries
        numread = fread(ja, sizeof(int), nz, ptr_myfile);
        if (numread != nz){
            fprintf(stderr," error in fread");
            return 1;
        }
        for (i = 0; i < nz; ++i){
            ja[i] = bswap_32(ja[i]);//swap between big endian and small endian
            //ja[i] = htobe32(ja[i]);//to big endian

            // printf("ja[i] = %d\n", ja[i]);//%f %p %s %c
        }

        double *a = malloc(nz*sizeof(*a)); //values of all non-zero entries
        numread = fread(a, sizeof(double), nz, ptr_myfile);
        if (numread != nz){
            fprintf(stderr," error in fread");
            return 1;
        }

        for (i = 0; i < nz; ++i){
            //        a[i] = bswap_64(a[i]);
            a[i] = bswap_double(a[i]);//need to be optimized
            //        double bigEndian_double(double a)
            //        printf("a[i] = %20.16e\n", a[i]);//%f %p %s %c

        }

        csptr csmat = malloc(sizeof(*csmat));
        setupCS(csmat,m,1);//int setupCS(csptr amat, int len, int job)
        int firstpos = 0;
        int *nnzrow = csmat->nnzrow;

        for (i = 0; i < n; ++i)
        {
            //ia-1;ja-1
            //printf("\n i=%d",i); getchar();
            nnzrow[i] = nnzptr[i];//get the lenth of each row
            csmat->pj[i] = (int *) Malloc(nnzrow[i]*sizeof(int), "colunms2csptr" );//ja=mat->ja[i]=vbmat->ja[i];
            memcpy(csmat->pj[i], &ja[firstpos], nnzrow[i]*sizeof(int));
            //memcpy(mat->pj[i],&ja[firstpos],nnzrow[i]*sizeof(int));
            csmat->pa[i] = (double*) Malloc(nnzrow[i]*sizeof(double), "colunms2csptr" );//ja=mat->ja[i]=vbmat->ja[i];
            memcpy(csmat->pa[i],&a[firstpos],nnzrow[i]*sizeof(double));
            //printf("\n vbmat->ja[i]=%d",vbmat->ja[i]); getchar();
            firstpos += nnzrow[i];
        }
        outputcsmat(csmat,"csmat.coo",1);//getchar();

    }

//    if(header == 1211214){
//        long int m;//length of the vector
//        numread = fread(&m, sizeof(int), 1, ptr_myfile);
//        m = bswap_32(m);//row
//        printf("m is %ld \n", m );


//        double *a = malloc(m*sizeof(*a)); //values of all non-zero entries
//        numread = fread(a, sizeof(double), m, ptr_myfile);
//        if (numread != m){
//            fprintf(stderr," error in fread");
//            return 1;
//        }

//        for (i = 0; i < m; ++i){
//            //        a[i] = bswap_64(a[i]);
//            a[i] = bswap_double(a[i]);//need to be optimized
//            //        double bigEndian_double(double a)
//            //        printf("a[i] = %20.16e\n", a[i]);//%f %p %s %c

//        }
//        output_dblvector("rhs.coo",a,0, m);

//    }

    return 0;
}


