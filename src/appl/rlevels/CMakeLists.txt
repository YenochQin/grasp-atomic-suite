add_executable(rlevels rlevels.f90)
target_link_libraries_Fortran(rlevels PRIVATE core)
install(TARGETS rlevels DESTINATION bin/ COMPONENT rlevels)
