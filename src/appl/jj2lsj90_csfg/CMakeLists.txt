# Source list follows the GRASPG Makefile (e157d9f).
add_executable(jj2lsj_csfg
    getmixblock.f90
    getmixblock_I.f90
    idigit.f90
    idigit_I.f90
    lval.f90
    lval_I.f90
    packLS.f90
    packLS_I.f90
    jj2lsj_data_1_C.f90
    jj2lsj_data_2_C.f90
    jj2lsj_data_3_C.f90
    jj2lsj_code.f90
    jj2lsj2K.f90
)
setup_fortran_modules(jj2lsj_csfg)
target_include_directories(jj2lsj_csfg PRIVATE "${Fortran_MODULE_DIRECTORY_root}/jj2lsj_csfg")
target_link_libraries_Fortran(jj2lsj_csfg PRIVATE 9290_csfg mod_csfg rang90 9290 mod)
install(TARGETS jj2lsj_csfg DESTINATION bin/)
