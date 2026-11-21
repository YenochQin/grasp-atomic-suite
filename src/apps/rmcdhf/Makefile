EXE=${GRASP}/bin/rmcdhf
LIBS=-L ${GRASP}/lib/ -ldavidson_eigensolver -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/davidson_eigensolver -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	serial_mpi_context.o \
	csf_density_matrix_element.o \
	 \
	direct_angular_coefficient.o \
	 \
	exchange_angular_coefficient.o \
	 \
	assemble_orbital_potential_coefficients.o \
	 \
	load_csf_basis_data.o \
	 \
	configure_debug_output.o \
	 \
	assemble_hamiltonian_block.o \
	 \
	build_exchange_radial_potential.o \
	 \
	build_direct_radial_potential.o \
	 \
	apply_off_diagonal_integrals.o \
	 \
	check_orbital_damping.o \
	 \
	calculate_lagrange_multipliers.o \
	 \
	connect_mcp_files.o \
	 \
	open_mixing_output.o \
	 \
	open_summary_file.o \
	 \
	build_radial_inhomogeneous_terms.o \
	 \
	build_variation_source_term.o \
	 \
	clear_radial_inhomogeneous_terms.o \
	 \
	apply_lagrange_constraints.o \
	 \
	build_orbital_potential.o \
	 \
	evaluate_radial_self_consistency.o \
	 \
	write_csf_weights.o \
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
	read_al_average_options.o \
	 \
	print_scf_results.o \
	 \
	read_eol_level_weights.o \
	 \
	read_eol_average_options.o \
	 \
	read_scf_options.o \
	 \
	write_hamiltonian_records.o \
	 \
	integrate_radial_outward.o \
	 \
	integrate_radial_inward.o \
	 \
	print_radial_wavefunction.o \
	 \
	basis_csf_parity.o \
	 \
	basis_csf_two_j_plus_one.o \
	 \
	diagonalize_hamiltonian_block.o \
	 \
	solve_hamiltonian_blocks.o \
	 \
	find_least_consistent_orbital.o \
	 \
	update_scf_mixing.o \
	 \
	energy_outside_bounds.o \
	 \
	update_trial_orbital_energy.o \
	 \
	write_scf_radial_wavefunctions.o \
	 \
	orthogonalize_to_more_consistent_orbitals.o \
	 \
	orthogonalize_same_kappa_orbitals.o \
	 \
	solve_radial_orbital.o \
	 \
	improve_orbital.o \
	 \
	run_scf_iterations.o \
	 \
	write_input_summary.o \
	 \
	rmcdhf_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
