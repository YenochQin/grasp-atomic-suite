add_executable(gj90
    compute_gj_only_I.f90
    setout_gj_I.f90
    ../rhfs90/gethfd_I.f90
    ../rhfs90/getmixblock_I.f90
    ../rhfs90/matelt_I.f90
    ../rhfs90/rinthf_I.f90
    ../rhfs90/setdbg_I.f90
    compute_gj_only.f90
    setout_gj.f90
    ../rhfs90/gethfd.f90
    ../rhfs90/getmixblock.f90
    ../rhfs90/matelt.f90
    ../rhfs90/rinthf.f90
    ../rhfs90/setdbg.f90
    gj90.f90
)
target_link_libraries_Fortran(gj90 PUBLIC rang90 mcp90 9290 mod)
install(TARGETS gj90 DESTINATION bin/)

if(TARGET mpiu90)
    add_executable(gj90_mpi
        compute_gj_only_mpi_I.f90
        setout_gj_I.f90
        ../rhfs90/gethfd_I.f90
        ../rhfs90/getmixblock_I.f90
        ../rhfs90/matelt_I.f90
        ../rhfs90/rinthf_I.f90
        ../rhfs90/setdbg_I.f90
        compute_gj_only_mpi.f90
        setout_gj.f90
        ../rhfs90/gethfd.f90
        ../rhfs90/getmixblock.f90
        ../rhfs90/matelt.f90
        ../rhfs90/rinthf.f90
        ../rhfs90/setdbg.f90
        gj90_mpi.f90
    )
    target_link_libraries_Fortran(gj90_mpi PUBLIC mpiu90 rang90 mcp90 9290 mod)
    target_include_directories(gj90_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(gj90_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(gj90_mpi PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
    )
    install(TARGETS gj90_mpi DESTINATION bin/)
endif()
