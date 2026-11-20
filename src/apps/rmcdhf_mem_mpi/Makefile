EXE=${GRASP}/bin/rmcdhf_mem_mpi
LIBS=-L ${GRASP}/lib/ -ldavidson_eigensolver -lgrasp_mpi -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/davidson_eigensolver -I ${GRASP}/src/lib/mpi -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	rmcdhf_mcp_memory.o \
	read_mcp_structure_from_memory.o \
	 \
	read_mcp_coefficients_from_memory.o \
	 \
	load_mcp_records_into_memory.o \
	 \
	read_mcp_structure_into_memory.o \
	read_mcp_coefficients_into_memory.o \
	query_total_system_memory.o \
	 \
	query_available_system_memory.o \
	 \
	cache_mcp_files.o \
	 \
	csf_density_matrix_element.o \
	 \
	direct_angular_coefficient.o \
	 \
	exchange_angular_coefficient.o \
	 \
	assemble_orbital_potential_coefficients.o \
	 \
	build_exchange_radial_potential.o \
	 \
	build_direct_radial_potential.o \
	 \
	apply_lagrange_constraints.o \
	 \
	apply_off_diagonal_integrals.o \
	 \
	build_orbital_potential_mpi.o \
	 \
	evaluate_radial_self_consistency.o \
	 \
	write_csf_weights.o \
	 \
	check_orbital_damping.o \
	 \
	damp_orbital_update.o \
	 \
	calculate_deferred_corrections.o \
	 \
	solve_radial_eigenvalue.o \
	 \
	write_scf_energies.o \
	 \
	finish_summary.o \
	 \
	estimate_orbital_energy.o \
	 \
	read_al_level_weights.o \
	 \
	read_al_average_options_mpi.o \
	 \
	print_scf_results.o \
	 \
	read_eol_level_weights.o \
	 \
	read_eol_average_options_mpi.o \
	 \
	read_scf_options_mpi.o \
	 \
	write_hamiltonian_records_mpi.o \
	 \
	build_radial_inhomogeneous_terms.o \
	 \
	build_variation_source_term.o \
	 \
	clear_radial_inhomogeneous_terms.o \
	 \
	integrate_radial_outward.o \
	 \
	integrate_radial_inward.o \
	 \
	print_radial_wavefunction.o \
	 \
	energy_outside_bounds.o \
	 \
	update_trial_orbital_energy.o \
	 \
	solve_radial_orbital.o \
	 \
	calculate_lagrange_multipliers_mpi.o \
	 \
	orthogonalize_to_more_consistent_orbitals.o \
	 \
	orthogonalize_same_kappa_orbitals.o \
	 \
	improve_orbital_mpi.o \
	 \
	basis_csf_parity.o \
	 \
	basis_csf_two_j_plus_one.o \
	 \
	lodcslmpiGG.o \
	lodcslmpiGG_I.o \
	lodcsh2GG.o \
	lodcsh2GG_I.o \
	find_least_consistent_orbital.o \
	 \
	update_scf_mixing_mpi.o \
	 \
	write_scf_radial_wavefunctions.o \
	 \
	load_csf_basis_data_mpi.o \
	 \
	configure_debug_output.o \
	 \
	configure_debug_output_mpi.o \
	 \
	assemble_hamiltonian_block.o \
	 \
	connect_mcp_files.o \
	 \
	open_mixing_output.o \
	 \
	open_summary_file.o \
	 \
	write_input_summary.o \
	 \
	diagonalize_hamiltonian_block_mpi.o \
	 \
	solve_hamiltonian_blocks_mpi.o \
	 \
	run_scf_iterations_mpi.o \
	 \
	rmcdhf_mpi_main.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
