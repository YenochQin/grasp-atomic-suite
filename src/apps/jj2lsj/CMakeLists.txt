# Ordinary CSF jj -> LS conversion, independent of jj2lsj_csfg.
option(GRASP_JJ2LSJ_OPENMP "Enable OpenMP LS-row conversion in jj2lsj_omp" ON)
if(GRASP_JJ2LSJ_OPENMP)
    find_package(OpenMP COMPONENTS Fortran)
endif()

# Persistent controller state is explicit; OpenMP supplies reentrant locals to
# the concurrent build. Do not restore global implicit SAVE compatibility flags.

# Recompile only the concurrent call closure, from the shared source files.
# Direct objects resolve these symbols before the legacy static archives. No
# second grasp_shared/libcore library or duplicate global module state is introduced.
set(jj2lsj_kernels)
foreach(name IN ITEMS
        wigner_9j wigner_9j_with_zero wigner_6j
        wigner_6j_with_j_one wigner_6j_with_j_two wigner_6j_with_j_three
        wigner_6j_with_j_three_halves wigner_6j_with_j_four wigner_6j_with_j_one_half
        racah_coefficients angular_triangle_and_parity wigner_6j_triangle_selection)
    list(APPEND jj2lsj_kernels "../../lib/spin_angular/${name}.f90")
endforeach()
list(APPEND jj2lsj_kernels "../../lib/core/csf_quantum_queries.f90")
add_library(jj2lsj_thread_kernels OBJECT ${jj2lsj_kernels})
target_link_libraries_Fortran(jj2lsj_thread_kernels PRIVATE spin_angular core grasp_shared)

add_library(jj2lsj_conversion STATIC
    jj2lsj_profile.f90 jj2lsj_postprocessing.f90
    read_mixing_blocks.f90
    ascii_digit_value.f90  orbital_letter_to_l.f90  pack_ls_subshell_label.f90
    ls_jj_coefficients_low_j.f90 ls_jj_coefficients_f6.f90 ls_jj_coefficients_f7.f90
    jj2lsj_conversion.f90)
setup_fortran_modules(jj2lsj_conversion)
target_include_directories(jj2lsj_conversion PRIVATE
    "${Fortran_MODULE_DIRECTORY_root}/jj2lsj_conversion")
target_link_libraries_Fortran(jj2lsj_conversion PUBLIC core spin_angular grasp_shared)
target_sources(jj2lsj_conversion INTERFACE $<TARGET_OBJECTS:jj2lsj_thread_kernels>)
if(GRASP_JJ2LSJ_OPENMP AND OpenMP_Fortran_FOUND)
    set_property(TARGET jj2lsj_conversion PROPERTY JJ2LSJ_OPENMP TRUE)
    target_link_libraries(jj2lsj_conversion PUBLIC OpenMP::OpenMP_Fortran)
    target_link_libraries(jj2lsj_thread_kernels PRIVATE OpenMP::OpenMP_Fortran)
elseif(GRASP_JJ2LSJ_OPENMP)
    message(WARNING "OpenMP Fortran unavailable: jj2lsj_omp uses its serial backend")
endif()

add_executable(jj2lsj_omp jj2lsj_main.f90)
target_link_libraries_Fortran(jj2lsj_omp PRIVATE jj2lsj_conversion)
install(TARGETS jj2lsj_omp DESTINATION bin/ COMPONENT jj2lsj_omp)
