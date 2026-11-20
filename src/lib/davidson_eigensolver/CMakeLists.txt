add_library(davidson_eigensolver STATIC
    davidson_projected_matrix.f90
    davidson_iteration.f90
    davidson_initial_subspace.f90
    davidson_orthonormalization.f90
    davidson_basis_transform.f90
    davidson_residual_expansion.f90
    davidson_subspace_restart.f90
    davidson_convergence.f90
    davidson_solver_control.f90
    davidson_eigensolver.f90
    davidson_sparse_initial_estimates.f90
)
setup_fortran_modules(davidson_eigensolver)
target_link_libraries_Fortran(davidson_eigensolver PRIVATE core grasp_shared)
install(TARGETS davidson_eigensolver DESTINATION lib/)
