add_executable(rhfs
    engouth.f90
    engouth_I.f90
    gethfd.f90
    gethfd_I.f90
    getmixblock.f90
    getmixblock_I.f90
    opt6_C.f90
    matelt.f90
    matelt_I.f90
    rinthf.f90
    rinthf_I.f90
    setdbg.f90
    setdbg_I.f90
    setsum.f90
    setsum_I.f90
    strsum.f90
    strsum_I.f90
    hfsgg.f90
    hfsgg_I.f90
    hfs92.f90
)
target_link_libraries_Fortran(rhfs PUBLIC rang90 mcp90 9290 mod)
install(TARGETS rhfs DESTINATION bin/)

if(TARGET mpiu90)
    add_executable(rhfs_mpi
        engouth.f90
        engouth_I.f90
        gethfd.f90
        gethfd_I.f90
        getmixblock.f90
        getmixblock_I.f90
        opt6_C.f90
        matelt.f90
        matelt_I.f90
        rinthf.f90
        rinthf_I.f90
        setdbg.f90
        setdbg_I.f90
        setsum.f90
        setsum_I.f90
        strsum.f90
        strsum_I.f90
        hfsgg_mpi.f90
        hfsgg_mpi_I.f90
        hfs92_mpi.f90
    )
    target_link_libraries_Fortran(rhfs_mpi PUBLIC mpiu90 rang90 mcp90 9290 mod)
    target_include_directories(rhfs_mpi PRIVATE ${MPI_Fortran_INCLUDE_PATH})
    target_link_libraries(rhfs_mpi PRIVATE ${MPI_Fortran_LIBRARIES})
    set_target_properties(rhfs_mpi PROPERTIES
      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
    )
    install(TARGETS rhfs_mpi DESTINATION bin/)
endif()
