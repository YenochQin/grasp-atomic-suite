EXE=${GRASP}/bin/rtransition_mpi
LIBS=-L ${GRASP}/lib/ -lbreit_angular_support -ldavidson_eigensolver -lrang90 -lgrasp_mpi -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/breit_angular_support -I ${GRASP}/src/lib/davidson_eigensolver -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/mpi -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	resize_transition_label_arrays.o \
	 \
	resize_transition_coefficient_arrays.o \
	 \
	resize_transition_matrix_arrays.o \
	 \
	resize_transition_integral_arrays.o \
	 \
	open_transition_angular_cache.o \
	 \
	build_transition_bessel_functions.o \
	 \
	calculate_biorthogonal_overlap.o \
	 \
	map_merged_orbitals.o \
	 \
	prepare_transition_mixing_coefficients.o \
	 \
	single_particle_transition_element.o \
	 \
	write_transition_csf_weights.o \
	 \
	write_transition_level_energies.o \
	 \
	read_state_filenames.o \
	 \
	read_multipole_options.o \
	 \
	read_final_radial_wavefunctions.o \
	 \
	read_initial_radial_wavefunctions.o \
	 \
	read_transition_options.o \
	 \
	right_csf_subshell_occupation.o \
	 \
	right_csf_parity.o \
	 \
	right_csf_two_j_plus_one.o \
	 \
	right_csf_intermediate_coupling.o \
	 \
	right_csf_subshell_quantum_number.o \
	 \
	load_initial_transition_basis.o \
	 \
	load_final_transition_basis.o \
	 \
	load_initial_lsj_labels.o \
	 \
	load_final_lsj_labels.o \
	 \
	sort_transition_angular_coefficients.o \
	 \
	read_transition_angular_coefficients.o \
	 \
	write_transition_angular_coefficients.o \
	 \
	merge_transition_basis_records.o \
	 \
	merge_transition_csf_bases.o \
	 \
	read_transition_mixing_coefficients.o \
	 \
	write_lsj_transition_probability.o \
	 \
	write_transition_probability.o \
	 \
	calculate_oscillator_strengths.o \
	 \
	write_input_summary.o \
	 \
	check_transition_mixing_coefficients.o \
	 \
	rtransition_mpi_main.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
