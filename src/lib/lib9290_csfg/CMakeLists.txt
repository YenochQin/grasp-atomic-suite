# Source list follows the GRASPG Makefile (e157d9f).
add_library(9290_csfg STATIC
    indexsym.f90
    indexsym_I.f90
    analabv.f90
    analabv_I.f90
    csfg_ispar.f90
    csfg_ispar_I.f90
    labpack.f90
    labpack_I.f90
    csfg_itjpo.f90
    csfg_itjpo_I.f90
    csfg_jcup.f90
    csfg_jcup_I.f90
    csfg_setqna.f90
    csfg_setqna_I.f90
    fictious_csf.f90
    findtype.f90
    printfictcsf.f90
    printlabelv.f90
    setlaborb.f90
    setlaborb_I.f90
    testfiles.f90
    testfiles_I.f90
    lodstate_csfg.f90
    lodstate_csfg_I.f90
    lodcsh_csfg.f90
    lodcsh_csfg_I.f90
    lodcsh2_csfg.f90
    lodcsh2_csfg_I.f90
    setcsll_csfg.f90
    setcsll_csfg_I.f90
    setcsh_csfg.f90
    setcsh_csfg_I.f90
    rcsfsymexpand.f90
    rcsfsymexpand_I.f90
)
setup_fortran_modules(9290_csfg)
target_include_directories(9290_csfg PRIVATE "${Fortran_MODULE_DIRECTORY_root}/9290_csfg")
target_link_libraries_Fortran(9290_csfg PRIVATE mod_csfg 9290 mod)
if(TARGET mpiu90)
    target_sources(9290_csfg PRIVATE lodcslmpi_csfg.f90 lodcslmpi_csfg_I.f90 cslhmpi_csfg.f90 cslhmpi_csfg_I.f90)
    target_link_libraries_Fortran(9290_csfg PRIVATE mpiu90)
endif()
install(TARGETS 9290_csfg DESTINATION lib/)
