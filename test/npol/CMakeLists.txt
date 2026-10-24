add_library(npol_test_support STATIC support.f90)
setup_fortran_modules(npol_test_support)
target_link_libraries_Fortran(npol_test_support PUBLIC mod)

function(npol_checked target)
    if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
        target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace)
    endif()
endfunction()

function(npol_test name)
    add_test(NAME ${name} COMMAND ${ARGN})
    set_tests_properties(${name} PROPERTIES LABELS npol_code TIMEOUT 60)
endfunction()

add_executable(npol_unit_driver unit_driver.f90)
target_link_libraries_Fortran(npol_unit_driver PRIVATE npol_test_support mod 9290)
npol_checked(npol_unit_driver)
foreach(mode kernel compose bounds)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${mode}")
    npol_test(npol_code_${mode} $<TARGET_FILE:npol_unit_driver> ${mode})
    set_tests_properties(npol_code_${mode} PROPERTIES WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${mode}")
endforeach()

if(NOT Python3_Interpreter_FOUND)
    message(WARNING "NP input/output/entrypoint and MPI wrappers require Python3; these tests are NOT registered.")
    return()
endif()
set(check "${CMAKE_CURRENT_SOURCE_DIR}/check.py")
foreach(mode input output)
    npol_test(npol_code_${mode} "${Python3_EXECUTABLE}" "${check}" ${mode} $<TARGET_FILE:npol_unit_driver>)
endforeach()
npol_test(npol_code_entrypoints "${Python3_EXECUTABLE}" "${check}" entrypoints
    "${PROJECT_SOURCE_DIR}" $<TARGET_FILE:npol_unit_driver>)
npol_test(npol_code_upstream_copy "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/check_upstream.py" "${PROJECT_SOURCE_DIR}")

foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(variant MATCHES "_mpi$" AND NOT TARGET mpiu90)
        continue()
    endif()
    set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
    set(target npol_consumers_${variant})
    add_executable(${target} consumers.F90 "${app}/ypot.f90" "${app}/ypot_I.f90"
        "${app}/dacon.f90" "${app}/dacon_I.f90" "${PROJECT_SOURCE_DIR}/src/lib/lib9290/rinti.f90")
    setup_fortran_modules(${target})
    target_link_libraries_Fortran(${target} PRIVATE npol_test_support mod 9290)
    npol_checked(${target})
    if(variant MATCHES "_mpi$")
        target_compile_definitions(${target} PRIVATE NPOL_MPI)
        target_link_libraries_Fortran(${target} PRIVATE mpiu90)
        target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
        if(MPIEXEC_EXECUTABLE)
            npol_test(npol_code_consumers_${variant} "${Python3_EXECUTABLE}" "${check}" run
                "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(npol_code_consumers_${variant} PROPERTIES PROCESSORS 2)
        endif()
    else()
        npol_test(npol_code_consumers_${variant} "${Python3_EXECUTABLE}" "${check}" run $<TARGET_FILE:${target}>)
    endif()
endforeach()

if(TARGET mpiu90 AND MPIEXEC_EXECUTABLE)
    add_executable(npol_mpi_driver mpi_driver.f90)
    target_link_libraries_Fortran(npol_mpi_driver PRIVATE npol_test_support mpiu90 mod 9290)
    target_link_libraries(npol_mpi_driver PRIVATE MPI::MPI_Fortran)
    npol_checked(npol_mpi_driver)
    foreach(ranks 1 2)
        npol_test(npol_code_mpi_broadcast_${ranks} "${Python3_EXECUTABLE}" "${check}" mpi
            "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:npol_mpi_driver> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(npol_code_mpi_broadcast_${ranks} PROPERTIES PROCESSORS ${ranks} TIMEOUT 120)
    endforeach()
    npol_test(npol_code_mpi_error_2 "${Python3_EXECUTABLE}" "${check}" mpi_error
        "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}
        $<TARGET_FILE:npol_mpi_driver> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(npol_code_mpi_error_2 PROPERTIES PROCESSORS 2 TIMEOUT 180)
endif()
