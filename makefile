# normally this files  need not be edited -- architecture-dependent
# make commands are in makefile.in -- 
# this makefile will make the library -- for testing PARARMS -- go to
# examples/grid or examples/general and see README files..
include makefile.in

# path of the header files of pARMS
ICFLAGS			=	-I./include

# path of the header files for implmentaion of pARMS
ISRCINC	                = -I./src/include	

# path of the header files
IFFLAGS			=  -I./include

# library path and directory declaration
LIB	=    ./lib/libparms.a
DIRS    =    ./ ./include ./lib ./src ./src/include ./src/DDPQ \
             ./src/FORTRAN ./examples ./examples/general \
             ./examples/grid ./examples/petsc ./examples/matrices  

# files from the src directory
OBJ1	=    ./src/parms_comm.o  \
	     ./src/parms_map.o  ./src/parms_mat.o  \
	     ./src/parms_mat_vcsr.o ./src/parms_mat_dvcsr.o  \
             ./src/parms_mem.o ./src/parms_operator.o  \
             ./src/parms_ilu_vcsr.o ./src/parms_pc.o \
	     ./src/parms_pc_bj.o  ./src/parms_pc_ras.o  \
	     ./src/parms_pc_schur.o  \
             ./src/parms_qsplit.o ./src/parms_solver.o  \
             ./src/parms_table.o ./src/parms_timer.o  \
             ./src/parms_vec.o ./src/parms_viewer.o ./src/fgmres.o \
             ./src/gmres.o ./src/parms_complex.o ./src/parms_vbilu_vcsr.o

# files from DDPQ directory
OBJ2 	=    ./src/DDPQ/arms2.o ./src/DDPQ/ilutpC.o  \
	     ./src/DDPQ/MatOps.o ./src/DDPQ/misc.o  \
	     ./src/DDPQ/PQ.o ./src/DDPQ/piluNEW.o  \
             ./src/DDPQ/setblks.o ./src/DDPQ/sets.o  \
	     ./src/DDPQ/svdInvC.o ./src/DDPQ/vbiluk.o \
	     ./src/DDPQ/vbilut.o ./src/DDPQ/barms2.o  \
	     ./src/DDPQ/barmsmisc.o ./src/DDPQ/vbiluNEW.o \
	     ./src/DDPQ/barmsmiscold.o ./src/DDPQ/cliques.o ./src/DDPQ/profile.o
#./src/DDPQ/ilut.o

# FORTRAN interface files
OBJ3	=   ./src/FORTRAN/parms_mapf.o ./src/FORTRAN/parms_matf.o \
            ./src/FORTRAN/parms_pcf.o ./src/FORTRAN/parms_solverf.o \
            ./src/FORTRAN/parms_timerf.o ./src/FORTRAN/parms_vecf.o \
            ./src/FORTRAN/parms_viewerf.o ./src/DDPQ/tools.o

OBJ     =     $(OBJ1) $(OBJ2) $(OBJ3)

default: $(LIB) 

all: $(LIB) tests

$(LIB): $(OBJ)
	if [ ! -d lib ]; then \
	  mkdir lib; \
	fi
	$(AR) $(ARFLAGS) $(LIB) $(OBJ) 

tests:  $(LIB)
	(cd examples/grid;make dd-grid.ex ) 
	(cd examples/general;make dd-HB-dse.ex)

petsc:  $(LIB)
	(cd examples/petsc;make dd-petsc.ex)
	(cd examples/petsc;make test.ex)


.c.o:
	${CC} ${ICFLAGS} ${ISRCINC} ${XIFLAGS} $(COPTFLAGS) \
	${CFLAGS} ${CFFLAGS} $< -c -o $@

.F.o:
	${F90} -FI ${IFFLAGS} ${FFLAGS} $< -c -o $(@F)

.f.o:
	${FC} ${FFLAGS} $< -c -o $@
#.c.o:
#	$(CC) $(CFLAGS) $(MPI_INC) $< -c -o $@

#.f.o:
#	$(FC) $(FFLAGS) $(MPI_INC) $< -c -o $@

cleanall: 
	@for dir in $(DIRS) ;\
          do \
          echo cleaning $$dir ;\
          (cd $$dir;  rm -rf *.a *.o *.ex* *core* out* \#*  *~) ;\
          done

cleanobj:
	@for dir in $(DIRS) ;\
	do \
          echo cleaning $$dir ;\
          (cd $$dir;  rm -rf  *.o  *core* out* sol.* \#*  *~) ;\
          done
