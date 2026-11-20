EXE=${GRASP}/bin/rangular_mpi
LIBS=-L ${GRASP}/lib/ -lrang90 -lbreit_angular_support -lgrasp_mpi -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/breit_angular_support -I ${GRASP}/src/lib/mpi -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	fndbeg.o \
	fndbeg_I.o \
	getinf.o \
	getinf_I.o \
	outsdampi.o \
	outsdampi_I.o \
	setdbg.o \
	setdbg_I.o \
	setdbgmpi.o \
	setdbgmpi_I.o \
	setmcp.o \
	setmcp_I.o \
	setmcpmpi.o \
	setmcpmpi_I.o \
	setsda.o \
	setsda_I.o \
	setsum.o \
	setsum_I.o \
	settmp.o \
	settmp_I.o \
	sort.o \
	sort_I.o \
	strsum.o \
	strsum_I.o \
	mcpmpi_gg.o \
	mcpmpi_gg_I.o \
	genmcpmpi.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
