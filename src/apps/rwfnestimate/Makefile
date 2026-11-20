EXE=${GRASP}/bin/rwfnestimate
LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	integrate_radial_tail.o \
	integrate_radial_tail_interface.o \
	integrate_homogeneous_radial_segment.o \
	integrate_homogeneous_radial_segment_interface.o \
	solve_homogeneous_radial_equation.o \
	solve_homogeneous_radial_equation_interface.o \
	estimate_hydrogenic_orbitals.o \
	estimate_hydrogenic_orbitals_interface.o \
	estimate_from_saved_orbitals.o \
	estimate_from_saved_orbitals_interface.o \
	estimate_thomas_fermi_orbitals.o \
	estimate_thomas_fermi_orbitals_interface.o \
	build_thomas_fermi_potential.o \
	build_thomas_fermi_potential_interface.o \
	print_uninitialized_orbitals.o \
	print_uninitialized_orbitals_interface.o \
	set_hydrogenic_screening.o \
	set_hydrogenic_screening_interface.o \
	configure_debug_output.o \
	configure_debug_output_interface.o \
	open_summary_file.o \
	open_summary_file_interface.o \
	write_input_summary.o \
	write_input_summary_interface.o \
	write_initial_wavefunction_summary.o \
	write_initial_wavefunction_summary_interface.o \
	write_initial_radial_wavefunctions.o \
	write_initial_radial_wavefunctions_interface.o \
	generate_initial_radial_wavefunctions.o \
	generate_initial_radial_wavefunctions_interface.o \
	read_calculation_options.o \
	read_calculation_options_interface.o \
	rwfnestimate_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
