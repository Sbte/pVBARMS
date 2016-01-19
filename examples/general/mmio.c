    /*
    *   Matrix Market I/O library for ANSI C
    *
    *   See http://math.nist.gov/MatrixMarket for details.
    *
    *
    */


    #include <stdio.h>
    #include <string.h>
    #include <stdlib.h>
    #include <ctype.h>
    #include <mpi.h>

    #include "mmio.h"


    int mm_is_valid(MM_typecode matcode)
    {
        if (!mm_is_matrix(matcode)) return 0;
        if (mm_is_dense(matcode) && mm_is_pattern(matcode)) return 0;
        if (mm_is_real(matcode) && mm_is_hermitian(matcode)) return 0;
        if (mm_is_pattern(matcode) && (mm_is_hermitian(matcode) ||
                                       mm_is_skew(matcode))) return 0;
        return 1;
    }

    int mm_read_banner(FILE *f, MM_typecode *matcode)
    {
        char line[MM_MAX_LINE_LENGTH];
        char banner[MM_MAX_TOKEN_LENGTH];
        char mtx[MM_MAX_TOKEN_LENGTH];
        char crd[MM_MAX_TOKEN_LENGTH];
        char data_type[MM_MAX_TOKEN_LENGTH];
        char storage_scheme[MM_MAX_TOKEN_LENGTH];
        char *p;


        mm_clear_typecode(matcode);

        if (fgets(line, MM_MAX_LINE_LENGTH, f) == NULL)
            return MM_PREMATURE_EOF;

        if (sscanf(line, "%s %s %s %s %s", banner, mtx, crd, data_type,
                   storage_scheme) != 5)
            return MM_PREMATURE_EOF;

        for (p=mtx; *p!='\0'; *p=tolower(*p),p++);  /* convert to lower case */
        for (p=crd; *p!='\0'; *p=tolower(*p),p++);
        for (p=data_type; *p!='\0'; *p=tolower(*p),p++);
        for (p=storage_scheme; *p!='\0'; *p=tolower(*p),p++);

        /* check for banner */
        if (strncmp(banner, MatrixMarketBanner, strlen(MatrixMarketBanner)) != 0)
            return MM_NO_HEADER;

        /* first field should be "mtx" */
        if (strcmp(mtx, MM_MTX_STR) != 0)
            return  MM_UNSUPPORTED_TYPE;
        mm_set_matrix(matcode);


        /* second field describes whether this is a sparse matrix (in coordinate
                storgae) or a dense array */


        if (strcmp(crd, MM_SPARSE_STR) == 0)
            mm_set_sparse(matcode);
        else
            if (strcmp(crd, MM_DENSE_STR) == 0)
                mm_set_dense(matcode);
            else
                return MM_UNSUPPORTED_TYPE;


        /* third field */

        if (strcmp(data_type, MM_REAL_STR) == 0)
            mm_set_real(matcode);
        else
            if (strcmp(data_type, MM_COMPLEX_STR) == 0)
                mm_set_complex(matcode);
            else
                if (strcmp(data_type, MM_PATTERN_STR) == 0)
                    mm_set_pattern(matcode);
                else
                    if (strcmp(data_type, MM_INT_STR) == 0)
                        mm_set_integer(matcode);
                    else
                        return MM_UNSUPPORTED_TYPE;


        /* fourth field */

        if (strcmp(storage_scheme, MM_GENERAL_STR) == 0)
            mm_set_general(matcode);
        else
            if (strcmp(storage_scheme, MM_SYMM_STR) == 0)
                mm_set_symmetric(matcode);
            else
                if (strcmp(storage_scheme, MM_HERM_STR) == 0)
                    mm_set_hermitian(matcode);
                else
                    if (strcmp(storage_scheme, MM_SKEW_STR) == 0)
                        mm_set_skew(matcode);
                    else
                        return MM_UNSUPPORTED_TYPE;


        return 0;
    }

    int mm_write_mtx_crd_size(FILE *f, int M, int N, int nz)
    {
        if (fprintf(f, "%d %d %d\n", M, N, nz) != 3)
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }

    int mm_read_mtx_crd_size(FILE *f, int *M, int *N, int *nz, MM_typecode matcode)
    {
        char line[MM_MAX_LINE_LENGTH];
        int num_items_read;

        /* set return null parameter values, in case we exit with errors */
        *M = *N = *nz = 0;

        /* now continue scanning until you reach the end-of-comments */
        do
        {
            if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
                return MM_PREMATURE_EOF;
        }while (line[0] == '%');

        /* line[] is either blank or has M,N, nz */
        if (sscanf(line, "%d %d %d", M, N, nz) != 3)
        {
            do
            {
                num_items_read = fscanf(f, "%d %d %d", M, N, nz);
                if (num_items_read == EOF) return MM_PREMATURE_EOF;
            }
            while (num_items_read != 3);
        }

        if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
            *nz = *nz * 2 - *N;
        else if (mm_is_skew(matcode))
            *nz = *nz * 2;

        return 0;
    }


    int mm_read_mtx_array_size(FILE *f, int *M, int *N)
    {
        char line[MM_MAX_LINE_LENGTH];
        int num_items_read;
        /* set return null parameter values, in case we exit with errors */
        *M = *N = 0;

        /* now continue scanning until you reach the end-of-comments */
        do
        {
            if (fgets(line,MM_MAX_LINE_LENGTH,f) == NULL)
                return MM_PREMATURE_EOF;
        }while (line[0] == '%');

        /* line[] is either blank or has M,N, nz */
        if (sscanf(line, "%d %d", M, N) == 2)
            return 0;

        else /* we have a blank line */
            do
        {
            num_items_read = fscanf(f, "%d %d", M, N);
            if (num_items_read == EOF) return MM_PREMATURE_EOF;
        }
        while (num_items_read != 2);

        return 0;
    }

    int mm_write_mtx_array_size(FILE *f, int M, int N)
    {
        if (fprintf(f, "%d %d\n", M, N) != 2)
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }

    int mm_read_array_data(FILE *f, int n, FLOAT *val, MM_typecode matcode)
    {
        int i;

        if (mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)
            double realvalue, imagvalue;

            for (i=0; i<n; i++){
                //            if (fscanf(f, "%lg %lg", &val[2*i], &val[2*i+1]) != 2)
                //                return MM_PREMATURE_EOF;
                if (fscanf(f, "%lg %lg", &realvalue, &imagvalue) != 2)
                    return MM_PREMATURE_EOF;
                val[i] = realvalue + imagvalue * I;
            }
    #endif
        }
        else if (mm_is_real(matcode))
        {
    #if defined(DBL)
            for (i=0; i<n; i++)
            {
                if (fscanf(f, "%lg\n", &val[i]) != 1)
                    return MM_PREMATURE_EOF;
            }
    #endif
        }
        else
            return MM_UNSUPPORTED_TYPE;
        return 0;

    }


    /*-------------------------------------------------------------------------*/

    /******************************************************************/
    /* use when I[], J[], and val[]J, and val[] are already allocated */
    /******************************************************************/

    int mm_read_mtx_crd_data(FILE *f, int M, int N, int nz, int ii[], int jj[],
                             FLOAT *val, MM_typecode matcode)
    {
        int i, j, fnz;
        fnz = nz;

        if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
            fnz = (fnz + N) / 2;
        else if (mm_is_skew(matcode))
            fnz = fnz / 2;

        if (mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)
            double realvalue, imagvalue;

            for (i=0; i<fnz; i++){
                //            if (fscanf(f, "%d %d %lg %lg", &ii[i], &jj[i], &val[2*i], &val[2*i+1])
                //                    != 4) return MM_PREMATURE_EOF;
                if (fscanf(f, "%d %d %lg %lg", &ii[i], &jj[i], &realvalue, &imagvalue) != 4)
                    return MM_PREMATURE_EOF;
                val[i] = realvalue + imagvalue * I;
            }
    #endif
        }
        else if (mm_is_real(matcode) || mm_is_integer(matcode))
        {
    #if defined(DBL)

            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg\n", &ii[i], &jj[i], &val[i]) != 3)
                    return MM_PREMATURE_EOF;

            }
    #endif
        }

        else if (mm_is_pattern(matcode))
        {
            for (i=0; i<fnz; i++){
                if (fscanf(f, "%d %d", &ii[i], &jj[i]) != 2)
                    return MM_PREMATURE_EOF;
                val[i] = 1.0;
            }
        }
        else
            return MM_UNSUPPORTED_TYPE;

        j = fnz;
        if ((mm_is_symmetric(matcode) && mm_is_real(matcode)) || (mm_is_hermitian(matcode) && mm_is_real(matcode)))
        {
            for (i=0; i<fnz; i++)
            {
                if (ii[i] == jj[i])
                    continue;
                ii[j] = jj[i];
                jj[j] = ii[i];
                val[j] = val[i];
                j += 1;
            }
        }
        else if (mm_is_hermitian(matcode) && mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)
            for (i=0; i<fnz; i++)
            {
                //            printf("i = %d, j = %d\n", i, j);//%f %p %s %c

                if (ii[i] == jj[i])
                    continue;
                ii[j] = jj[i];
                jj[j] = ii[i];
                val[j] = conj(val[i]);
                //            val[j*2] = val[i*2];
                //            val[j*2+1] = -val[i*2+1];
                j += 1;
            }
    #endif
        }
        else if (mm_is_symmetric(matcode) && mm_is_complex(matcode))
        {
            for (i=0; i<fnz; i++)
            {
                if (ii[i] == jj[i])
                    continue;
                ii[j] = jj[i];
                jj[j] = ii[i];
                val[j] = val[i];
                //            val[j*2] = val[i*2];
                //            val[j*2+1] = val[i*2+1];
                j += 1;
            }
        }
        else if (mm_is_skew(matcode))
        {
            for (i=0; i<fnz; i++)
            {
                ii[j] = jj[i];
                jj[j] = ii[i];
                val[j] = - val[i];
                j += 1;
            }
        }

        return 0;

    }

    int mm_partial_read_mtx_crd_data(FILE *f, int M, int N, int nz, int *ii, int *jj,
                                     FLOAT *val, int *idom, int *dom, int *perm, int *nB, int
                                     nBlock, MM_typecode matcode)
    {
        int i, fnz;
        int pid;

        fnz = nz;

        MPI_Comm_rank(MPI_COMM_WORLD, &pid);
        int *p = (int *) calloc(M+1, sizeof(int));
        int *bsz = (int *) calloc(nBlock+1, sizeof(int));

        int i1, i3;
        int i2 = 0;

        for (i1 = 0; i1 < nBlock; ++i1)
            bsz[i1+1] = bsz[i1] + nB[i1];
        int istart = idom[pid];
        int iend = idom[pid+1];

        for (i1 = istart; i1 < iend; ++i1)
        {
            i2 = dom[i1-1]-1; //index of permuted thing that is on this proc
            for (i3 = bsz[i2]; i3 < bsz[i2+1]; ++i3)
                p[i3] = 1;
        }

        if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
            fnz = (fnz + N) / 2;
        else if (mm_is_skew(matcode))
            fnz = fnz / 2;

        i3 = 0;
        int ir, jr;
        double vr1, vr2;

        // TODO: Doesn't work for symmetric/hermitian non-real matrices
        if (mm_is_complex(matcode))
        {
            printf("matrix is complex.\n");//%f %p %s %c

    #if defined(DBL_CMPLX)
            //        double realvalue, imagvalue;

            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg %lg", &ir, &jr, &vr1, &vr2)
                        != 4) return MM_PREMATURE_EOF;
                //            if (fscanf(f, "%d %d %lg %lg", &ii[i], &jj[i], &realvalue, &imagvalue) != 4)
                //                return MM_PREMATURE_EOF;
                //            val[i] = realvalue + imagvalue * I;

                //            if (p[perm[ir-1]])
                //            {
                //                ii[i3] = ir;
                //                jj[i3] = jr;
                //                val[2*i3] = vr1;
                //                val[2*i3+1] = vr2;
                //                i3++;
                //            }

                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1 + vr2 * I;
                    i3++;
                }
            }
    #endif

        }
        else if ((mm_is_symmetric(matcode) && mm_is_real(matcode)) || (mm_is_hermitian(matcode) && mm_is_real(matcode)))
        {
    //        printf(" in symm mode\n");//%f %p %s %c
    #if defined(DBL)

            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg\n", &ir, &jr, &vr1)
                        != 3)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1;
                    i3++;
                }
                if (p[perm[jr-1]] && ir != jr)
                {
                    ii[i3] = jr;
                    jj[i3] = ir;
                    val[i3] = vr1;
                    i3++;
                }
            }
    #endif

        }
        else if (mm_is_real(matcode))
        {
    //        printf(" in real mode");//%f %p %s %c

    #if defined(DBL)
            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg\n", &ir, &jr, &vr1)
                        != 3) return MM_PREMATURE_EOF;


                ////            if((val[i] >= -EPSILON) && (val[i] <= EPSILON))
                //            if((0 >= -EPSILON) && (0 <= EPSILON))
                ////            if(val[i] == 0.0)
                //            {
                //                printf("i is %d, j is %d.\n", ir, jr);//%f %p %s %c

                //            }
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1;
                    i3++;
                }

            }
    #endif
        }
        else if (mm_is_pattern(matcode))
        {
            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d", &ir, &jr)
                        != 2) return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    i3++;
                }
            }
        }
        else
            return MM_UNSUPPORTED_TYPE;

        printf("i3 %d fnz %d n %d\n", i3, fnz, M);

        free(p);
        free(bsz);

        return 0;

    }


    int mm_partial_read_mtx_crd_data_new(FILE *f, int M, int N, int nz, int *ii, int *jj,
                                         FLOAT *val, int *idom, int *dom, int *perm, int *nB, int
                                         nBlock, MM_typecode matcode)
    {
        int i, fnz;
        int pid;

        fnz = nz;

        MPI_Comm_rank(MPI_COMM_WORLD, &pid);
        int *p = (int *) calloc(M+1, sizeof(int));
        int *bsz = (int *) calloc(nBlock+1, sizeof(int));

        int i1, i3;
        int i2 = 0;

        for (i1 = 0; i1 < nBlock; ++i1)
            bsz[i1+1] = bsz[i1] + nB[i1];
        int istart = idom[pid];
        int iend = idom[pid+1];

        for (i1 = istart; i1 < iend; ++i1)
        {
            i2 = dom[i1-1]-1; //index of permuted thing that is on this proc
            for (i3 = bsz[i2]; i3 < bsz[i2+1]; ++i3)
                p[i3] = 1;
        }

        if (mm_is_symmetric(matcode) || mm_is_hermitian(matcode))
            fnz = (fnz + N) / 2;
        else if (mm_is_skew(matcode))
            fnz = fnz / 2;

        i3 = 0;
        int ir, jr;
        double vr1, vr2;

        if ((mm_is_symmetric(matcode) && mm_is_real(matcode)) || (mm_is_hermitian(matcode) && mm_is_real(matcode)))
        {
    //        printf(" in symm mode\n");//%f %p %s %c
    #if defined(DBL)

            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg\n", &ir, &jr, &vr1)
                        != 3)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1;
                    i3++;
                }
                if (p[perm[jr-1]] && ir != jr)
                {
                    ii[i3] = jr;
                    jj[i3] = ir;
                    val[i3] = vr1;
                    i3++;
                }
            }
    #endif

        }
        else if (mm_is_hermitian(matcode) && mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)

            for (i=0; i<fnz; i++)
            {
    //            if (fscanf(f, "%d %d %lg %lg", &ii[i], &jj[i], &realvalue, &imagvalue) != 4)

                if (fscanf(f, "%d %d %lg %lg\n", &ir, &jr, &vr1, &vr2) != 4)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1 + vr2 * I;
                    i3++;
                }
                if (p[perm[jr-1]] && ir != jr)
                {
                    ii[i3] = jr;
                    jj[i3] = ir;
                    val[i3] = vr1 - vr2 * I;
                    i3++;
                }
            }
    #endif
        }
        else if (mm_is_symmetric(matcode) && mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)

            for (i=0; i<fnz; i++)
            {

                if (fscanf(f, "%d %d %lg %lg\n", &ir, &jr, &vr1, &vr2) != 4)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1 + vr2 * I;
                    i3++;
                }
                if (p[perm[jr-1]] && ir != jr)
                {
                    ii[i3] = jr;
                    jj[i3] = ir;
                    val[i3] =  vr1 + vr2 * I;
                    i3++;
                }
            }
    #endif
        }
        else if (mm_is_skew(matcode) && mm_is_real(matcode))
        {
    #if defined(DBL)

            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg\n", &ir, &jr, &vr1)
                        != 3)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1;
                    i3++;
                }
                if (p[perm[jr-1]] && ir != jr)
                {
                    ii[i3] = jr;
                    jj[i3] = ir;
                    val[i3] = -vr1;
                    i3++;
                }
            }
    #endif
        }
        else if (mm_is_skew(matcode) && mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)

            for (i=0; i<fnz; i++)
            {

                if (fscanf(f, "%d %d %lg %lg\n", &ir, &jr, &vr1, &vr2) != 4)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1 + vr2 * I;
                    i3++;
                }
                if (p[perm[jr-1]] && ir != jr)
                {
                    ii[i3] = jr;
                    jj[i3] = ir;
                    val[i3] = - vr1 - vr2 * I;
                    i3++;
                }
            }
    #endif
        }
        else if (mm_is_real(matcode))
        {
    //        printf(" in real mode");//%f %p %s %c

    #if defined(DBL)
            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d %lg\n", &ir, &jr, &vr1)
                        != 3) return MM_PREMATURE_EOF;


                ////            if((val[i] >= -EPSILON) && (val[i] <= EPSILON))
                //            if((0 >= -EPSILON) && (0 <= EPSILON))
                ////            if(val[i] == 0.0)
                //            {
                //                printf("i is %d, j is %d.\n", ir, jr);//%f %p %s %c

                //            }
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1;
                    i3++;
                }

            }
    #endif
        }
        else if (mm_is_complex(matcode))
        {
    #if defined(DBL_CMPLX)

            for (i=0; i<fnz; i++)
            {

                if (fscanf(f, "%d %d %lg %lg\n", &ir, &jr, &vr1, &vr2) != 4)
                    return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    val[i3] = vr1 + vr2 * I;
                    i3++;
                }
            }
    #endif
        }
        else if (mm_is_pattern(matcode))
        {
            for (i=0; i<fnz; i++)
            {
                if (fscanf(f, "%d %d", &ir, &jr)
                        != 2) return MM_PREMATURE_EOF;
                if (p[perm[ir-1]])
                {
                    ii[i3] = ir;
                    jj[i3] = jr;
                    i3++;
                }
            }
        }
        else
            return MM_UNSUPPORTED_TYPE;

    //    printf("i3 %d fnz %d n %d\n", i3, fnz, M);

        free(p);
        free(bsz);

        return 0;

    }
    //int mm_read_mtx_crd_entry(FILE *f, int *I, int *J,
    //                          double *real, double *imag, MM_typecode matcode)
    //{
    //    if (mm_is_complex(matcode))
    //    {
    //        if (fscanf(f, "%d %d %lg %lg", I, J, real, imag)
    //                != 4) return MM_PREMATURE_EOF;
    //    }
    //    else if (mm_is_real(matcode))
    //    {
    //        if (fscanf(f, "%d %d %lg\n", I, J, real)
    //                != 3) return MM_PREMATURE_EOF;

    //    }

    //    else if (mm_is_pattern(matcode))
    //    {
    //        if (fscanf(f, "%d %d", I, J) != 2) return MM_PREMATURE_EOF;
    //    }
    //    else
    //        return MM_UNSUPPORTED_TYPE;

    //    return 0;

    //}


    /************************************************************************
        mm_read_mtx_crd()  fills M, N, nz, array of values, and return
                            type code, e.g. 'MCRS'

                            if matrix is complex, values[] is of size 2*nz,
                                (nz pairs of real/imaginary values)
    ************************************************************************/

    int mm_read_mtx_crd(char *fname, int *M, int *N, int *nz, int **ii, int **jj,
                        FLOAT **val, MM_typecode *matcode)
    {
        int ret_code;
        FILE *f;

        if (strcmp(fname, "stdin") == 0) f=stdin;
        else
            if ((f = fopen(fname, "r")) == NULL)
                return MM_COULD_NOT_READ_FILE;


        if ((ret_code = mm_read_banner(f, matcode)) != 0)
            return ret_code;

        if (!(mm_is_valid(*matcode) && mm_is_sparse(*matcode) &&
              mm_is_matrix(*matcode)))
            return MM_UNSUPPORTED_TYPE;

        if ((ret_code = mm_read_mtx_crd_size(f, M, N, nz, *matcode)) != 0)
            return ret_code;


        *ii = (int *)  malloc(*nz * sizeof(int));
        *jj = (int *)  malloc(*nz * sizeof(int));
        *val = NULL;

        if (mm_is_complex(*matcode))
        {
            *val = (FLOAT *) malloc(*nz * 2 * sizeof(FLOAT));
            ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *ii, *jj, *val,
                                            *matcode);
            if (ret_code != 0) return ret_code;
        }
        else if (mm_is_real(*matcode))
        {
            *val = (FLOAT *) malloc(*nz * sizeof(FLOAT));
            ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *ii, *jj, *val,
                                            *matcode);
            if (ret_code != 0) return ret_code;
        }

        else if (mm_is_pattern(*matcode))
        {
            ret_code = mm_read_mtx_crd_data(f, *M, *N, *nz, *ii, *jj, *val,
                                            *matcode);
            if (ret_code != 0) return ret_code;
        }

        if (f != stdin) fclose(f);
        return 0;
    }

    int mm_write_banner(FILE *f, MM_typecode matcode)
    {
        char *str = mm_typecode_to_str(matcode);
        int ret_code;

        ret_code = fprintf(f, "%s %s\n", MatrixMarketBanner, str);
        free(str);
        if (ret_code !=2 )
            return MM_COULD_NOT_WRITE_FILE;
        else
            return 0;
    }


    /**
    *  Create a new copy of a string s.  mm_strdup() is a common routine, but
    *  not part of ANSI C, so it is included here.  Used by mm_typecode_to_str().
    *
    */
    char *mm_strdup(const char *s)
    {
        int len = strlen(s);
        char *s2 = (char *) malloc((len+1)*sizeof(char));
        return strcpy(s2, s);
    }

    char  *mm_typecode_to_str(MM_typecode matcode)
    {
        char buffer[MM_MAX_LINE_LENGTH];
        char *types[4];
        char *mm_strdup(const char *);
        //    int error = 0;

        /* check for MTX type */
        //    if (mm_is_matrix(matcode))
        types[0] = MM_MTX_STR;
        //    else
        //        error=1;

        /* check for CRD or ARR matrix */
        if (mm_is_sparse(matcode))
            types[1] = MM_SPARSE_STR;
        else
            if (mm_is_dense(matcode))
                types[1] = MM_DENSE_STR;
            else
                return NULL;

        /* check for element data type */
        if (mm_is_real(matcode))
            types[2] = MM_REAL_STR;
        else
            if (mm_is_complex(matcode))
                types[2] = MM_COMPLEX_STR;
            else
                if (mm_is_pattern(matcode))
                    types[2] = MM_PATTERN_STR;
                else
                    if (mm_is_integer(matcode))
                        types[2] = MM_INT_STR;
                    else
                        return NULL;


        /* check for symmetry type */
        if (mm_is_general(matcode))
            types[3] = MM_GENERAL_STR;
        else
            if (mm_is_symmetric(matcode))
                types[3] = MM_SYMM_STR;
            else
                if (mm_is_hermitian(matcode))
                    types[3] = MM_HERM_STR;
                else
                    if (mm_is_skew(matcode))
                        types[3] = MM_SKEW_STR;
                    else
                        return NULL;

        sprintf(buffer,"%s %s %s %s", types[0], types[1], types[2], types[3]);
        return mm_strdup(buffer);

    }


    //            char Tau_string[80];

    //            sprintf(Tau_string, "%0.2f", prm->eps);

    //            printf("Tau_string value is %s\n", Tau_string);//%f %p %s %c

    //            char str1[80];
    //            strcpy(str1, "perm_");

    //            char* perm_file_name = strcat(str1, curname);
    //            char str2[80];
    //            strcpy(str2, "nB_");
    //            char* nB_file_name = strcat(str2, curname);

    //            printf("perm_file_name value is %s\n", perm_file_name);//%f %p %s %c

    //            perm_file_name = strcat(perm_file_name, "_");
    //            perm_file_name = strcat(perm_file_name, Tau_string);
    //            perm_file_name = strcat(perm_file_name, ".txt");

    //            nB_file_name = strcat(nB_file_name, "_");
    //            nB_file_name = strcat(nB_file_name, Tau_string);
    //            nB_file_name = strcat(nB_file_name, ".txt");

    //            printf("curname value is %s\n", curname);//%f %p %s %c

    //            output_intvector(perm_file_name, perm, 0, csmat->n);
    //            output_intvector(nB_file_name, nB, 0, nBlock);


    //            goto nextmat;
