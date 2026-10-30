# Source list follows the GRASPG Makefile (e157d9f).
add_library(mod_csfg STATIC
    csfg_parameters_M.f90
    csfg_decide_C.f90
    symexpand_mod.f90
    symmatrix_restart_C.f90
    symmatrix_mod.f90
    csfg_memory_man.f90
    csfg_tv_C.f90
    csfg_scf_C.f90
)
setup_fortran_modules(mod_csfg)
target_include_directories(mod_csfg PRIVATE "${Fortran_MODULE_DIRECTORY_root}/mod_csfg")
target_link_libraries_Fortran(mod_csfg PRIVATE mod)
install(TARGETS mod_csfg DESTINATION lib/)
