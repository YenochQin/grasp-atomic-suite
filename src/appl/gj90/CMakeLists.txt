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
