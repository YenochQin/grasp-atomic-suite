# Source list follows the GRASPG Makefile (e157d9f).
add_executable(rmixaccumulate_csfg
    rmixaccumulate_csfg.f90
)
setup_fortran_modules(rmixaccumulate_csfg)
target_include_directories(rmixaccumulate_csfg PRIVATE "${Fortran_MODULE_DIRECTORY_root}/rmixaccumulate_csfg")
target_link_libraries_Fortran(rmixaccumulate_csfg PRIVATE 9290_csfg mod_csfg 9290 mod)
install(TARGETS rmixaccumulate_csfg DESTINATION bin/)
