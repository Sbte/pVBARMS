SHELL	= /bin/sh

include ../../makefile.in


# path of the header files of pARMS
ICFLAGS			=	-I../../include

# path of the header files for implmentaion of pARMS
ISRCINC	                = -I../../src/include	

# path of the header files
IFFLAGS			=  -I../../include

# pARMS lib and its directory
#LIBDIR			=	$../../lib
#LIB			=	${LIBDIR}/libparms.a
#LIBFLAGS		=	-L${LIBDIR}
#SHFLAGS			=	
PARMS_LIBS		=	-L../../lib -lparms


allexe: Bdd-HB-dse.ex Bdd-mtx-dse.ex Bdd-bin-dse.ex


Bdd-mtx-dse.ex: Bdd-mtx-dse.o aux.o skit.o mmio.o
	${LINKER} ${CLINKFLAGS} -o $@ Bdd-mtx-dse.o aux.o skit.o mmio.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR} ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB}

Bdd-bin-dse.ex: Bdd-bin-dse.o aux.o skit.o mmio.o
	${LINKER} ${CLINKFLAGS} -o $@ Bdd-bin-dse.o aux.o skit.o mmio.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR} ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB}


Bdd-HB-dse.ex: Bdd-HB-dse.o aux.o skit.o
	${LINKER} ${CLINKFLAGS} -o $@ Bdd-HB-dse.o aux.o skit.o \
	${SHFLAGS} ${LIBFLAGS} ${PARMS_LIBS} \
	${MPI_LIB_DIR} ${MPI_LIB} ${LAPACK_BLAS_LIB_DIR} ${LAPACK_BLAS_LIB} \
	${MATH_LIB_DIR} ${MATH_LIB} 

## the next test run uses metis - link to metis must be provided
# METIS header files 
XIFLAGS	= -I/home/saad/scicom/METIS-4.0/Lib
# METIS library path
XLIB	= -L/home/saad/scicom/METIS-4.0 -lmetis




.c.o:
	${CC} ${ICFLAGS} ${ISRCINC} ${XIFLAGS} $(COPTFLAGS) \
	${CFLAGS} ${CFFLAGS} $< -c -o $@

.F.o:
	${F90} ${IFFLAGS} ${FFLAGS} $< -c -o $(@F)

.f.o:
	${FC} ${FFLAGS} $< -c -o $(@F)

clean:
	${RM} -rf *.o core *~ *# *.ex paramfile.*




