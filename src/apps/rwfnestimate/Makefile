EXE=${GRASP}/bin/rwfnestimate
LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	integrate_radial_tail.o \
	 \
	integrate_homogeneous_radial_segment.o \
	 \
	solve_homogeneous_radial_equation.o \
	 \
	estimate_hydrogenic_orbitals.o \
	 \
	estimate_from_saved_orbitals.o \
	 \
	estimate_thomas_fermi_orbitals.o \
	 \
	build_thomas_fermi_potential.o \
	 \
	print_uninitialized_orbitals.o \
	 \
	set_hydrogenic_screening.o \
	 \
	configure_debug_output.o \
	 \
	open_summary_file.o \
	 \
	write_input_summary.o \
	 \
	write_initial_wavefunction_summary.o \
	 \
	write_initial_radial_wavefunctions.o \
	 \
	generate_initial_radial_wavefunctions.o \
	 \
	read_calculation_options.o \
	 \
	rwfnestimate_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
