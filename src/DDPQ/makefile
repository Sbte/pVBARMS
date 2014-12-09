# normally this files  need not be edited -- architecture-dependent
# make commands are in makefile.in -- 
# this makefile will make the library -- for testing PARARMS -- go to
# examples/grid or examples/general and see README files..
include ../../makefile.in

# path of the header files of pARMS
ICFLAGS			=	-I../../include

# path of the header files for implmentaion of pARMS
ISRCINC	                = -I../include	

# path of the header files
IFFLAGS			=  -I../../include

OBJ = cliques.o cliquestest.o sets.o

cliquestest: $(OBJ)
	${LINKER} ${CLINKFLAGS} -o $@ ${OBJ} \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR} ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} 

.c.o:
	${CC} -Wall ${ICFLAGS} ${ISRCINC} ${XIFLAGS} $(COPTFLAGS) \
	${CFLAGS} ${CFFLAGS} $< -c -o $@
