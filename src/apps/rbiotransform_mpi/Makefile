EXE=${GRASP}/bin/rbiotransform_mpi
LIBS=-L ${GRASP}/lib/ -lrang90 -lgrasp_mpi -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/mpi -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	open_biorthogonal_angular_cache.o \
	open_biorthogonal_angular_cache_interface.o \
	invert_banded_matrix.o \
	invert_banded_matrix_interface.o \
	calculate_orbital_overlap_matrix.o \
	calculate_orbital_overlap_matrix_interface.o \
	combine_vectors.o \
	combine_vectors_interface.o \
	write_matrix.o \
	write_matrix_interface.o \
	copy_vector.o \
	copy_vector_interface.o \
	read_state_filenames.o \
	read_state_filenames_interface.o \
	read_state_mixing_coefficients.o \
	read_state_mixing_coefficients_interface.o \
	calculate_subshell_overlap.o \
	calculate_subshell_overlap_interface.o \
	sum_integer_vector.o \
	sum_integer_vector_interface.o \
	find_integer_extremum.o \
	find_integer_extremum_interface.o \
	vector_inner_product.o \
	vector_inner_product_interface.o \
	invert_dense_submatrix.o \
	invert_dense_submatrix_interface.o \
	final_state_radial_integral.o \
	final_state_radial_integral_interface.o \
	initial_state_radial_integral.o \
	initial_state_radial_integral_interface.o \
	interpolate_final_orbital.o \
	interpolate_final_orbital_interface.o \
	interpolate_initial_orbital.o \
	interpolate_initial_orbital_interface.o \
	prepare_biorthogonal_subshells.o \
	prepare_biorthogonal_subshells_interface.o \
	lodcslBio.o \
	lodcslBio_I.o \
	read_final_radial_wavefunctions.o \
	read_final_radial_wavefunctions_interface.o \
	read_initial_radial_wavefunctions.o \
	read_initial_radial_wavefunctions_interface.o \
	print_matrix_symmetry.o \
	print_matrix_symmetry_interface.o \
	factor_lower_upper.o \
	factor_lower_upper_interface.o \
	fill_vector.o \
	fill_vector_interface.o \
	multiply_submatrices.o \
	multiply_submatrices_interface.o \
	sort_biorthogonal_mcp.o \
	sort_biorthogonal_mcp_interface.o \
	read_biorthogonal_mcp.o \
	read_biorthogonal_mcp_interface.o \
	write_biorthogonal_mcp.o \
	write_biorthogonal_mcp_interface.o \
	build_subshell_transformation_matrix.o \
	build_subshell_transformation_matrix_interface.o \
	write_transformed_radial_wavefunctions.o \
	write_transformed_radial_wavefunctions_interface.o \
	configure_radial_grid.o \
	configure_radial_grid_interface.o \
	scale_vector.o \
	scale_vector_interface.o \
	load_biorthogonal_csf_basis_mpi.o \
	load_biorthogonal_csf_basis_mpi_interface.o \
	load_merged_csf_basis.o \
	load_merged_csf_basis_interface.o \
	copy_state_csf_data.o \
	copy_state_csf_data_interface.o \
	apply_subshell_transfer.o \
	apply_subshell_transfer_interface.o \
	apply_occupation_factor.o \
	apply_occupation_factor_interface.o \
	transpose_matrix.o \
	transpose_matrix_interface.o \
	combine_upper_lower_factors.o \
	combine_upper_lower_factors_interface.o \
	transform_ci_coefficients.o \
	transform_ci_coefficients_interface.o \
	transform_orbital_sets.o \
	transform_orbital_sets_interface.o \
	rbiotransform_mpi_main.o

$(EXE): $(OBJS)
	$(FC_MPI) -o $@ $? $(FC_MPILD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC_MPI) -c $(FC_MPIFLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
