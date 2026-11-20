EXE=${GRASP}/bin/rci
LIBS=-L ${GRASP}/lib/ -lrang90 -lbreit_angular_support -l9290 -ldavidson_eigensolver -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/breit_angular_support -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/davidson_eigensolver -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	write_ci_eigenvectors.o \
	radial_density_ratio.o \
	ci_restart_units.o \
	collect_breit_angular_coefficient.o \
	build_breit_bessel_kernels.o \
	 \
	multiply_sparse_hamiltonian.o reduce_sparse_matrix_product.o \
	 \
	write_level_energies.o \
	 \
	fullerton_rinker_k_function.o \
	 \
	fullerton_rinker_l_function.o \
	 \
	mohr_self_energy_factor.o \
	 \
	klarsfeld_maquet_self_energy_factor.o \
	 \
	hydrogenic_self_energy_factor.o \
	 \
	accumulate_closed_shell_breit.o \
	evaluate_breit_hartree_z.o \
	 \
	coulomb_radial_integral.o \
	 \
	lookup_slater_integral.o \
	 \
	transverse_radial_integral.o \
	 \
	breit_radial_integral.o \
	 \
	lookup_breit_radial_integral.o \
	 \
	breit_radial_integral_1.o \
	 \
	breit_radial_integral_2.o \
	 \
	breit_radial_integral_3.o \
	 \
	breit_radial_integral_4.o \
	 \
	breit_radial_integral_5.o \
	 \
	breit_radial_integral_6.o \
	 \
	satisfies_breit_group1_triangle.o \
	 \
	satisfies_breit_group2_triangle.o \
	 \
	tabulate_breit_integrals_1.o \
	 \
	tabulate_breit_integrals_2.o \
	 \
	tabulate_slater_integrals.o \
	 \
	read_ci_options.o \
	 \
	write_hamiltonian_records.o \
	 \
	hydrogenic_radial_overlap.o \
	 \
	one_electron_integral.o \
	 \
	find_breit_integral_index.o \
	 \
	initialize_dense_eigenvectors.o \
	 \
	initialize_sparse_eigenvectors.o \
	 \
	kinetic_energy_integral.o \
	 \
	read_ci_mixing_coefficients.o \
	 \
	read_ci_restart.o \
	 \
	multiply_restart_hamiltonian.o \
	 \
	diagonalize_hamiltonian_block.o \
	 \
	build_nuclear_charge_distribution.o \
	 \
	build_qed_self_energy.o \
	 \
	tabulate_self_energy_integrals.o \
	 \
	load_csf_basis_data.o \
	 \
	configure_ci_debug_output.o \
	 \
	open_mixing_output.o \
	 \
	connect_ci_restart_file.o \
	 \
	open_summary_file.o \
	 \
	electron_screening_factor.o \
	 \
	write_input_summary.o \
	 \
	satisfies_coulomb_triangle.o \
	 \
	build_second_order_vacuum_polarization.o \
	 \
	build_fourth_order_vacuum_polarization.o \
	 \
	build_vacuum_polarization.o \
	 \
	vinti_radial_integral.o \
	 \
	calculate_vinti_integral.o \
	 \
	vacuum_polarization_integral.o \
	 \
	calculate_vacuum_polarization_integral.o \
	read_ci_level_selection.o \
	 \
	assemble_hamiltonian_block.o \
	 \
	assemble_ci_matrix.o \
	 \
	assemble_ci_matrix_from_restart.o \
	 \
	prepare_correction_integral_tables.o \
	 \
	solve_hamiltonian_blocks.o \
	 \
	rci_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
