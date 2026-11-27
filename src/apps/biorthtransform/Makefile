EXE=${GRASP}/bin/biorthtransform
LIBS=-L ${GRASP}/lib/ -lrang90 -l9290 -lgrasp_shared
FC_MODULES= -I ${GRASP}/src/lib/spin_angular -I ${GRASP}/src/lib/core -I ${GRASP}/src/lib/grasp_shared

OBJS= \
	open_biorthogonal_angular_cache.o \
	 \
	invert_banded_matrix.o \
	 \
	calculate_orbital_overlap_matrix.o \
	 \
	write_matrix.o \
	 \
	scale_vector.o \
	 \
	copy_vector.o \
	 \
	select_biorthogonal_symmetries.o \
	 \
	read_state_mixing_coefficients.o \
	 \
	calculate_subshell_overlap.o \
	 \
	sum_integer_vector.o \
	 \
	find_integer_extremum.o \
	 \
	final_state_radial_integral.o \
	 \
	initial_state_radial_integral.o \
	 \
	vector_inner_product.o \
	 \
	interpolate_final_orbital.o \
	 \
	interpolate_initial_orbital.o \
	 \
	invert_dense_submatrix.o \
	 \
	prepare_biorthogonal_subshells.o \
	 \
	read_final_radial_wavefunctions.o \
	 \
	read_initial_radial_wavefunctions.o \
	 \
	print_matrix_symmetry.o \
	factor_lower_upper.o \
	 \
	fill_vector.o \
	 \
	multiply_submatrices.o \
	 \
	sort_biorthogonal_mcp.o \
	 \
	read_biorthogonal_mcp.o \
	 \
	write_biorthogonal_mcp.o \
	 \
	build_subshell_transformation_matrix.o \
	 \
	write_transformed_radial_wavefunctions.o \
	 \
	configure_radial_grid.o \
	 \
	load_merged_csf_basis.o \
	 \
	copy_state_csf_data.o \
	 \
	apply_subshell_transfer.o \
	 \
	apply_occupation_factor.o \
	 \
	transpose_matrix.o \
	 \
	combine_upper_lower_factors.o \
	 \
	combine_vectors.o \
	 \
	prepare_biorthogonal_mcp.o \
	 \
	transform_ci_coefficients.o \
	 \
	transform_orbital_sets.o \
	 \
	biorthtransform_main.o

$(EXE): $(OBJS)
	$(FC) -o $@ $? $(FC_LD) $(LIBS) $(LAPACK_LIBS)

%.o: %.f90
	$(FC) -c $(FC_FLAGS) $(FC_MODULES) -o $@ $<

clean:
	-@rm $(EXE)
	-@rm *.o *.mod
