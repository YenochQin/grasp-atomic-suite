if(GRASP_BUILD_SERIAL_APPS)
    add_executable(gj
        compute_gj_only_I.f90
        setout_gj_I.f90
        ../rhfs/gethfd_I.f90
        ../rhfs/getmixblock_I.f90
        ../rhfs/matelt_I.f90
        ../rhfs/rinthf_I.f90
        ../rhfs/setdbg_I.f90
        compute_gj_only.f90
        setout_gj.f90
        ../rhfs/gethfd.f90
        ../rhfs/getmixblock.f90
        ../rhfs/matelt.f90
        ../rhfs/rinthf.f90
        ../rhfs/setdbg.f90
        gj.f90
    )
    target_link_libraries_Fortran(gj PUBLIC rang90 mcp90 9290 mod)
    install(TARGETS gj DESTINATION bin/)
endif()

if(TARGET mpiu90)
    add_executable(gj_mpi
        compute_gj_only_mpi_I.f90
        setout_gj_I.f90
        ../rhfs/gethfd_I.f90
        ../rhfs/getmixblock_I.f90
        ../rhfs/matelt_I.f90
        ../rhfs/rinthf_I.f90
        ../rhfs/setdbg_I.f90
        compute_gj_only_mpi.f90
        setout_gj.f90
        ../rhfs/gethfd.f90
        ../rhfs/getmixblock.f90
        ../rhfs/matelt.f90
        ../rhfs/rinthf.f90
        ../rhfs/setdbg.f90
        gj_mpi.f90
    )
    target_link_libraries_Fortran(gj_mpi PUBLIC mpiu90 rang90 mcp90 9290 mod)
    target_include_directories(gj_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(gj_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(gj_mpi PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
    )
    install(TARGETS gj_mpi DESTINATION bin/)
endif()
