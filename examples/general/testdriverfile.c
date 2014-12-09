//#include <stdio.h>
////#include <process.h>
//#include <string.h>
//#include <stdlib.h>

//int main()
//{
//    FILE *stream;
////    char msg[]="this is a test";
//    int a[3] = {1, 2 ,3};
//    char buf[20];
//    if ((stream=fopen("dummy.fil","w+"))==NULL)
//    {
//        fprintf(stderr,"cannot open output file.\n");
//        return 1;
//    }
//    /*write some data to the file*/
////    fwrite(msg,1,strlen(msg)+1,stream);
//    fwrite(a,sizeof(int),3 ,stream);
//    /*seek to the beginning of the file*/
//    fseek(stream,0,SEEK_SET);
//    /*read the data and display it*/
////    fread(buf,1,strlen(msg)+1,stream);

//    int *b = malloc(3*sizeof(*b));
//    fread(b,sizeof(int),3,stream);
//    printf("b = %d, %d, %d\n", b[0], b[1], b[2]);
//    fclose(stream);
//    free(b);
////    system("pause");
//    return 0;
//}



//#include<endian.h>
#include<stdio.h>
//#include<iostream>
//#include<stdint.h>
#include <stdlib.h>
#include<byteswap.h>//biswap_32
//using namespace std;
typedef struct parms_vcsr {
  int      n;         	//!< the dimension of the matrix
  int      *nnzrow;    //!< the length of each row
  int      *space;     //!<  length of space ( a work array)
  int      off_proc_n; //!< size of off-processor contributions
  /*!
        Parameter: pj = An indirect pointer to store column indices.
   */
  int      **pj;
  /*! parameter: pa = An indirect pointer to store corresponding nonzero entries.
  */
  double    **pa;
} *parms_vcsr, *csptr;

//double double_swap(double value){
//    union v {
//        double       f;
//        uint64_t    i;
//    };
//        v val = value;
//    val.i = bswap_64(val.i);
//    return val.f;
//}
//uint64_t bigEndian_double(double a)
double bswap_double(double a)
{
//    uint64_t b;
    double b;
    unsigned char *src = (unsigned char *)&a,
                  *dst = (unsigned char *)&b;
//    if (is_littleEndian())
//    {
    int i;
    for (i = 0; i < 8; ++i) {
        dst[i] = src[7-i];
    }
//        dst[0] = src[7];
//        dst[1] = src[6];
//        dst[2] = src[5];
//        dst[3] = src[4];
//        dst[4] = src[3];
//        dst[5] = src[2];
//        dst[6] = src[1];
//        dst[7] = src[0];
//    }
//    else
//        b = *(uint64_t *)&a;
    return b;

}


int main()
{
    //    int counter;
    FILE *ptr_myfile;
    //    struct rec my_record;
    //    char buffer[20];
    long int header;
    int i;//for loop

    ptr_myfile = fopen("mat000.bin","r");//r and rb are the same
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



    long int m ;    // = header(1); //row;
    long int n ;  //   = header(2); //column;
    long int nz;//     = header(3);%nnz

//    nnz = fread(fd,m,'int32');  %nonzeros per row
//    sum_nz = sum(nnz);
//    long int m =
    fread(&header, sizeof(int), 1, ptr_myfile);
    header = bswap_32(header);
    printf("header is %ld \n", header );

    fread(&m, sizeof(int), 1, ptr_myfile);
    m = bswap_32(m);//row
    printf("m is %ld \n", m );

    fread(&n, sizeof(int), 1, ptr_myfile);
    n = bswap_32(n);//column
    printf("n is %ld \n", n );


    fread(&nz, sizeof(int), 1, ptr_myfile);
    nz = bswap_32(nz);// nnz
    printf("nz is %ld \n", nz);

    int *nnzptr = malloc(m*sizeof(*nnzptr)); //nnz of each row
//    int *nnz = (int*) malloc (m*sizeof(int));

    long int numread;
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
        printf("a[i] = %20.16e\n", a[i]);//%f %p %s %c

    }

    csptr csmat = malloc(sizeof(*csmat));
//    printf("a[0] = %20.16e\n", a[0]);//%f %p %s %c

//    long int sum_nz = 0;
//    for (i = 0; i < m; ++i){
//        nnzptr[i] = htobe32(nnzptr[i]);
//        sum_nz += nnzptr[i];
////        printf("nnz[i] = %d\n", nnz[i]);//%f %p %s %c
//    }

//    for (i = 0; i < m; ++i){
//        nnz[i] = htobe32(nnz[i]);
//        printf("nnz[i] = %d\n", nnz[i]);//%f %p %s %c
//    }
//    fread(nnz, sizeof(int), m, ptr_myfile);
//numread=fread(list,sizeof(char),25,stream);
//    for (i = 0; i < m; ++i)
//        printf("nnz[i] = %d\n", nnz[i]);//%f %p %s %c
    //uint32_t htobe32(uint32_t host_32bits);           // host to big-endian encoding

//    for (size_t idx = 0; idx != 3; ++idx) {
//        fread(&header, sizeof(int), 1, ptr_myfile);
//        header = htobe32(header);
//    }

//    long int *nnz = malloc(n*sizeof(*nnz));;//     = header(3);%nnz



    //cout <<"header is " << header << '\n';
//        printf("header is %ld\n", header);//%f %p %s %c
    //fread(&my_record,sizeof(struct rec),1,ptr_myfile);
    //fread (void *__restrict __ptr, size_t __size, size_t __n, FILE *__restrict __stream)
    //    fclose(ptr_myfile);
    return 0;
}
