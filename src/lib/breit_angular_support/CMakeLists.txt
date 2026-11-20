add_library(breit_angular_support STATIC
    # Keep the collector in a separate archive member so RCI can supply its own definition.
    breit_radial_integral_coefficients.f90
    breit_angular_coefficient_buffer.f90
    collect_breit_angular_coefficient.f90
)
setup_fortran_modules(breit_angular_support)
target_link_libraries_Fortran(breit_angular_support PRIVATE grasp_shared core)
install(TARGETS breit_angular_support DESTINATION lib/)
