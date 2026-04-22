EXE=${GRASP}/bin/gj90
LIBS=-L ${GRASP}/lib/ -lrang90 -lmcp90 -l9290 -lmod
FC_MODULES= -I ${GRASP}/src/lib/librang90 -I ${GRASP}/src/lib/libmcp90 -I ${GRASP}/src/lib/lib9290 -I ${GRASP}/src/lib/libmod
RHFS90_DIR=../rhfs90

OBJS= \
	compute_gj_only_I.o \
	setout_gj_I.o \
	gethfd_I.o \
	getmixblock_I.o \
	matelt_I.o \
	rinthf_I.o \
	setdbg_I.o \
	compute_gj_only.o \
	setout_gj.o \
	gethfd.o \
	getmixblock.o \
	matelt.o \
	rinthf.o \
	setdbg.o \
	gj90.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

$(filter gethfd_I.o getmixblock_I.o matelt_I.o rinthf_I.o setdbg_I.o gethfd.o getmixblock.o matelt.o rinthf.o setdbg.o,$(OBJS)): %.o: $(RHFS90_DIR)/%.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
