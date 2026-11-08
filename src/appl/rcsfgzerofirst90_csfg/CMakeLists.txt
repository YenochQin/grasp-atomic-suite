# Source list follows the GRASPG Makefile (e157d9f).
add_executable(rcsfgzerofirst_csfg
    symzf_mod.f90
    removeblank_I.f90
    lodcsl_Part_I.f90
    lodcsl_Zero_I.f90
    lodcsl_Full_I.f90
    set_CSF_ZFlist_I.f90
    set_CSF_number_I.f90
    lodcsl_Part.f90
    lodcsl_Zero.f90
    lodcsl_Full.f90
    set_CSF_ZFlist.f90
    set_CSF_number.f90
    removeblank.f90
    set_CSF_numberF_I.f90
    set_CSF_numberF.f90
    RCSFzerofirst.f90
)
setup_fortran_modules(rcsfgzerofirst_csfg)
target_include_directories(rcsfgzerofirst_csfg PRIVATE "${Fortran_MODULE_DIRECTORY_root}/rcsfgzerofirst_csfg")
target_link_libraries_Fortran(rcsfgzerofirst_csfg PRIVATE 9290_csfg mod_csfg 9290 mod)
install(TARGETS rcsfgzerofirst_csfg DESTINATION bin/)
