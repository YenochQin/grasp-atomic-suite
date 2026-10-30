# Keep the executable names used by the original GRASPG workflows.
foreach(program rcsfgexpand_csfg rcsfgextend_csfg rcsfgsplit_csfg
                rdistHmatrix_csfg rcsfgblocksplit_csfg rlevels_csfg)
    add_executable(${program} ${program}.f90)
    setup_fortran_modules(${program})
    target_link_libraries_Fortran(${program} PRIVATE mod_csfg mod)
    install(TARGETS ${program} DESTINATION bin/)
endforeach()
install(PROGRAMS rsave_csfg DESTINATION bin/)
