# Ordinary CSF jj -> LS conversion, independent of jj2lsj_csfg.
option(GRASP_JJ2LSJ_OPENMP "Enable OpenMP LS-row conversion in jj2lsj_omp" ON)
if(GRASP_JJ2LSJ_OPENMP)
    find_package(OpenMP COMPONENTS Fortran)
endif()

# Persistent controller state is explicit; OpenMP supplies reentrant locals to
# the concurrent build. Do not restore global implicit SAVE compatibility flags.

# Recompile only the concurrent call closure, from the shared source files.
# Direct objects resolve these symbols before the legacy static archives. No
# second libmod/lib9290 library or duplicate global module state is introduced.
set(jj2lsj_kernels)
foreach(name nine nine0 sixj sixj1 sixj2 sixj3 sixj35 sixj4 sixj5 dracah ittk ixjtik)
    list(APPEND jj2lsj_kernels "../../lib/librang90/${name}.f90")
endforeach()
foreach(name iq jqs itjpo ispar)
    list(APPEND jj2lsj_kernels "../../lib/lib9290/${name}.f90")
endforeach()
add_library(jj2lsj_thread_kernels OBJECT ${jj2lsj_kernels})
target_link_libraries_Fortran(jj2lsj_thread_kernels PRIVATE rang90 9290 mod)

add_library(jj2lsj_conversion STATIC
    jj2lsj_profile.f90 jj2lsj_postprocess.f90
    getmixblock.f90 getmixblock_I.f90
    idigit.f90 idigit_I.f90 lval.f90 lval_I.f90 packLS.f90 packLS_I.f90
    jj2lsj_data_1_C.f90 jj2lsj_data_2_C.f90 jj2lsj_data_3_C.f90
    jj2lsj_code.f90)
setup_fortran_modules(jj2lsj_conversion)
target_include_directories(jj2lsj_conversion PRIVATE
    "${Fortran_MODULE_DIRECTORY_root}/jj2lsj_conversion")
target_link_libraries_Fortran(jj2lsj_conversion PUBLIC 9290 rang90 mod)
target_sources(jj2lsj_conversion INTERFACE $<TARGET_OBJECTS:jj2lsj_thread_kernels>)
if(GRASP_JJ2LSJ_OPENMP AND OpenMP_Fortran_FOUND)
    set_property(TARGET jj2lsj_conversion PROPERTY JJ2LSJ_OPENMP TRUE)
    target_link_libraries(jj2lsj_conversion PUBLIC OpenMP::OpenMP_Fortran)
    target_link_libraries(jj2lsj_thread_kernels PRIVATE OpenMP::OpenMP_Fortran)
elseif(GRASP_JJ2LSJ_OPENMP)
    message(WARNING "OpenMP Fortran unavailable: jj2lsj_omp uses its serial backend")
endif()

add_executable(jj2lsj_omp jj2lsj2K.f90)
target_link_libraries_Fortran(jj2lsj_omp PRIVATE jj2lsj_conversion)
install(TARGETS jj2lsj_omp DESTINATION bin/ COMPONENT jj2lsj_omp)
