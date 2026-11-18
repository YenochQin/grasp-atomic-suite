add_executable(rwfnestimate
    rwfnestimate_input.f90
    tail.f90
    tail_I.f90
    sbstep.f90
    sbstep_I.f90
    solvh.f90
    solvh_I.f90
    frmhyd.f90
    frmhyd_I.f90
    frmrwf.f90
    frmrwf_I.f90
    frmtfp.f90
    frmtfp_I.f90
    tfpot.f90
    tfpot_I.f90
    prtrem.f90
    prtrem_I.f90
    screenpar.f90
    screenpar_I.f90
    setdbg.f90
    setdbg_I.f90
    setsum.f90
    setsum_I.f90
    strsum.f90
    strsum_I.f90
    summry.f90
    summry_I.f90
    wrtrwf.f90
    wrtrwf_I.f90
    genrwf.f90
    genrwf_I.f90
    getinfo.f90
    getinf_I.f90
    erwf.f90
)
target_link_libraries_Fortran(rwfnestimate PUBLIC rang90 9290 mod)

# Keep the header parser and discard the unused Parquet code in the Rust archive.
include(CheckLinkerFlag)
if(APPLE)
    set(rwfnestimate_dead_strip "LINKER:-dead_strip")
else()
    set(rwfnestimate_dead_strip "LINKER:--gc-sections")
endif()
check_linker_flag(Fortran "${rwfnestimate_dead_strip}" GRASP_RWFNESTIMATE_DEAD_STRIP)
if(GRASP_RWFNESTIMATE_DEAD_STRIP)
    target_link_options(rwfnestimate PRIVATE "${rwfnestimate_dead_strip}")
endif()

install(TARGETS rwfnestimate DESTINATION bin/)
