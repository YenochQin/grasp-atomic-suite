EXE=${GRASP}/bin/rtransition
LIBS=-L ${GRASP}/lib/ -ldavidson_eigensolver -lrang90 -lbreit_angular_support -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/davidson_eigensolver -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/breit_angular_support -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	resize_transition_label_arrays.o \
	resize_transition_label_arrays_interface.o \
	resize_transition_coefficient_arrays.o \
	resize_transition_coefficient_arrays_interface.o \
	resize_transition_matrix_arrays.o \
	resize_transition_matrix_arrays_interface.o \
	resize_transition_integral_arrays.o \
	resize_transition_integral_arrays_interface.o \
	open_transition_angular_cache.o \
	open_transition_angular_cache_interface.o \
	build_transition_bessel_functions.o \
	build_transition_bessel_functions_interface.o \
	calculate_biorthogonal_overlap.o \
	calculate_biorthogonal_overlap_interface.o \
	map_merged_orbitals.o \
	map_merged_orbitals_interface.o \
	prepare_transition_mixing_coefficients.o \
	prepare_transition_mixing_coefficients_interface.o \
	single_particle_transition_element.o \
	single_particle_transition_element_interface.o \
	write_transition_csf_weights.o \
	write_transition_csf_weights_interface.o \
	write_transition_level_energies.o \
	write_transition_level_energies_interface.o \
	read_state_filenames.o \
	read_state_filenames_interface.o \
	read_final_radial_wavefunctions.o \
	read_final_radial_wavefunctions_interface.o \
	read_initial_radial_wavefunctions.o \
	read_initial_radial_wavefunctions_interface.o \
	read_multipole_options.o \
	read_multipole_options_interface.o \
	read_transition_options.o \
	read_transition_options_interface.o \
	right_csf_subshell_occupation.o \
	right_csf_subshell_occupation_interface.o \
	right_csf_parity.o \
	right_csf_parity_interface.o \
	right_csf_two_j_plus_one.o \
	right_csf_two_j_plus_one_interface.o \
	right_csf_intermediate_coupling.o \
	right_csf_intermediate_coupling_interface.o \
	right_csf_subshell_quantum_number.o \
	right_csf_subshell_quantum_number_interface.o \
	load_initial_transition_basis.o \
	load_initial_transition_basis_interface.o \
	load_final_transition_basis.o \
	load_final_transition_basis_interface.o \
	load_initial_lsj_labels.o \
	load_initial_lsj_labels_interface.o \
	load_final_lsj_labels.o \
	load_final_lsj_labels_interface.o \
	lodcslm.o \
	lodcslm_I.o \
	read_transition_angular_coefficients.o \
	read_transition_angular_coefficients_interface.o \
	sort_transition_angular_coefficients.o \
	sort_transition_angular_coefficients_interface.o \
	write_transition_angular_coefficients.o \
	write_transition_angular_coefficients_interface.o \
	merge_transition_basis_records.o \
	merge_transition_basis_records_interface.o \
	merge_transition_csf_bases.o \
	merge_transition_csf_bases_interface.o \
	read_transition_mixing_coefficients.o \
	read_transition_mixing_coefficients_interface.o \
	write_lsj_transition_probability.o \
	write_lsj_transition_probability_interface.o \
	write_transition_probability.o \
	write_transition_probability_interface.o \
	calculate_oscillator_strengths.o \
	calculate_oscillator_strengths_interface.o \
	setcslm.o \
	setcslm_I.o \
	write_input_summary.o \
	write_input_summary_interface.o \
	check_transition_mixing_coefficients.o \
	check_transition_mixing_coefficients_interface.o \
	rtransition_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
