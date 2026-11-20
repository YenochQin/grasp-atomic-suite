EXE=${GRASP}/bin/rmcdhf_mpi
LIBS=-L ${GRASP}/lib/ -ldavidson_eigensolver -lgrasp_mpi -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/davidson_eigensolver -I ${GRASP}/src/lib/mpi -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	csf_density_matrix_element.o \
	csf_density_matrix_element_interface.o \
	direct_angular_coefficient.o \
	direct_angular_coefficient_interface.o \
	exchange_angular_coefficient.o \
	exchange_angular_coefficient_interface.o \
	assemble_orbital_potential_coefficients.o \
	assemble_orbital_potential_coefficients_interface.o \
	build_exchange_radial_potential.o \
	build_exchange_radial_potential_interface.o \
	build_direct_radial_potential.o \
	build_direct_radial_potential_interface.o \
	apply_lagrange_constraints.o \
	apply_lagrange_constraints_interface.o \
	apply_off_diagonal_integrals.o \
	apply_off_diagonal_integrals_interface.o \
	build_orbital_potential_mpi.o \
	build_orbital_potential_mpi_interface.o \
	evaluate_radial_self_consistency.o \
	evaluate_radial_self_consistency_interface.o \
	write_csf_weights.o \
	write_csf_weights_interface.o \
	check_orbital_damping.o \
	check_orbital_damping_interface.o \
	damp_orbital_update.o \
	damp_orbital_update_interface.o \
	calculate_deferred_corrections.o \
	calculate_deferred_corrections_interface.o \
	solve_radial_eigenvalue.o \
	solve_radial_eigenvalue_interface.o \
	write_scf_energies.o \
	write_scf_energies_interface.o \
	finish_summary.o \
	finish_summary_interface.o \
	estimate_orbital_energy.o \
	estimate_orbital_energy_interface.o \
	read_al_level_weights.o \
	read_al_level_weights_interface.o \
	read_al_average_options_mpi.o \
	read_al_average_options_mpi_interface.o \
	print_scf_results.o \
	print_scf_results_interface.o \
	read_eol_level_weights.o \
	read_eol_level_weights_interface.o \
	read_eol_average_options_mpi.o \
	read_eol_average_options_mpi_interface.o \
	read_scf_options_mpi.o \
	read_scf_options_mpi_interface.o \
	write_hamiltonian_records_mpi.o \
	write_hamiltonian_records_mpi_interface.o \
	build_radial_inhomogeneous_terms.o \
	build_radial_inhomogeneous_terms_interface.o \
	build_variation_source_term.o \
	build_variation_source_term_interface.o \
	clear_radial_inhomogeneous_terms.o \
	clear_radial_inhomogeneous_terms_interface.o \
	integrate_radial_outward.o \
	integrate_radial_outward_interface.o \
	integrate_radial_inward.o \
	integrate_radial_inward_interface.o \
	print_radial_wavefunction.o \
	print_radial_wavefunction_interface.o \
	energy_outside_bounds.o \
	energy_outside_bounds_interface.o \
	update_trial_orbital_energy.o \
	update_trial_orbital_energy_interface.o \
	solve_radial_orbital.o \
	solve_radial_orbital_interface.o \
	calculate_lagrange_multipliers_mpi.o \
	calculate_lagrange_multipliers_mpi_interface.o \
	orthogonalize_to_more_consistent_orbitals.o \
	orthogonalize_to_more_consistent_orbitals_interface.o \
	orthogonalize_same_kappa_orbitals.o \
	orthogonalize_same_kappa_orbitals_interface.o \
	improve_orbital_mpi.o \
	improve_orbital_mpi_interface.o \
	basis_csf_parity.o \
	basis_csf_parity_interface.o \
	basis_csf_two_j_plus_one.o \
	basis_csf_two_j_plus_one_interface.o \
	lodcslmpiGG.o \
	lodcslmpiGG_I.o \
	lodcsh2GG.o \
	lodcsh2GG_I.o \
	find_least_consistent_orbital.o \
	find_least_consistent_orbital_interface.o \
	update_scf_mixing_mpi.o \
	update_scf_mixing_mpi_interface.o \
	write_scf_radial_wavefunctions.o \
	write_scf_radial_wavefunctions_interface.o \
	load_csf_basis_data_mpi.o \
	load_csf_basis_data_mpi_interface.o \
	configure_debug_output.o \
	configure_debug_output_interface.o \
	configure_debug_output_mpi.o \
	configure_debug_output_mpi_interface.o \
	assemble_hamiltonian_block.o \
	assemble_hamiltonian_block_interface.o \
	connect_mcp_files.o \
	connect_mcp_files_interface.o \
	open_mixing_output.o \
	open_mixing_output_interface.o \
	open_summary_file.o \
	open_summary_file_interface.o \
	write_input_summary.o \
	write_input_summary_interface.o \
	diagonalize_hamiltonian_block_mpi.o \
	diagonalize_hamiltonian_block_mpi_interface.o \
	solve_hamiltonian_blocks_mpi.o \
	solve_hamiltonian_blocks_mpi_interface.o \
	run_scf_iterations_mpi.o \
	run_scf_iterations_mpi_interface.o \
	rmcdhf_mpi_main.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
