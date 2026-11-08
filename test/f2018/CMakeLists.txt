add_executable(test.f2018_davidson davidson_contracts.F90)
target_link_libraries_Fortran(test.f2018_davidson PRIVATE dvd90 9290 mod)
target_link_libraries(test.f2018_davidson PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME f2018_davidson_contracts COMMAND test.f2018_davidson)
add_executable(test.f2018_nmtejj nmtejj_probe.f90)
add_executable(test.f2018_linear_algebra linear_algebra_contracts.f90)
target_link_libraries_Fortran(test.f2018_linear_algebra PRIVATE mod)
target_link_libraries(test.f2018_linear_algebra PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME f2018_linear_algebra_contracts COMMAND test.f2018_linear_algebra)
add_executable(test.f2018_storage_bytes storage_bytes.f90)
target_link_libraries_Fortran(test.f2018_storage_bytes PRIVATE mod)
add_test(NAME f2018_storage_bytes COMMAND test.f2018_storage_bytes)
target_link_libraries_Fortran(test.f2018_nmtejj PRIVATE rang90 mod)
if(Python3_Interpreter_FOUND)
    add_test(NAME f2018_nmtejj_parity
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_nmtejj.py
        $<TARGET_FILE:test.f2018_nmtejj>)
endif()

if(Python3_Interpreter_FOUND)
    foreach(driver acceptance first_sweep)
        add_test(NAME f2018_${driver}_driver COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_${driver}.py)
    endforeach()
endif()

# Exercise the actual MATRIX entry points with controlled assembly/eigensolver
# seams. Each executable has its own production modules and saved FIRST state.
foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi rmcdhf_mem_strict_mpi rmcdhf_mem_diag_mpi)
    if(variant MATCHES "_mpi$" AND NOT TARGET mpiu90)
        continue()
    endif()
    set(directory ${variant})
    if(variant MATCHES "^rmcdhf_mem_(strict|diag)_mpi$")
        string(REPLACE "rmcdhf_" "rmcdhf90_" directory ${variant})
    endif()
    set(app "${PROJECT_SOURCE_DIR}/src/appl/${directory}")
    set(matrix matrix)
    set(eigen maneig)
    if(variant MATCHES "_mpi$")
        set(matrix matrixmpi)
        set(eigen maneigmpi)
    endif()
    set(target test.f2018_ci_${variant})
    add_executable(${target} ci_damping.F90 "${app}/${matrix}.f90" "${app}/${matrix}_I.f90"
        "${app}/${eigen}_I.f90" "${app}/setham_I.f90")
    set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
    target_link_libraries_Fortran(${target} PRIVATE mod)
    target_link_libraries(${target} PRIVATE ${BLAS_LIBRARIES})
    target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace)
    if(variant MATCHES "mem")
        target_compile_definitions(${target} PRIVATE CI_MEM)
        target_sources(${target} PRIVATE "${app}/rmcdhf_mem_C.f90" "${app}/read03_mem_I.f90" "${app}/read1_mem_I.f90")
    endif()
    if(variant MATCHES "_mpi$")
        target_compile_definitions(${target} PRIVATE CI_MPI)
        target_link_libraries_Fortran(${target} PRIVATE mpiu90)
        foreach(ranks 1 2 4)
            add_test(NAME f2018_ci_damping_${variant}_${ranks}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
            set_tests_properties(f2018_ci_damping_${variant}_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
        endforeach()
    else()
        target_sources(${target} PRIVATE "${app}/mpi_s.f90")
        add_test(NAME f2018_ci_damping_${variant} COMMAND ${target})
        set_tests_properties(f2018_ci_damping_${variant} PROPERTIES TIMEOUT 30)
    endif()
endforeach()

if(TARGET mpiu90 AND Python3_Interpreter_FOUND)
    add_executable(test.f2018_mpi_contracts mpi_contracts.f90)
    target_link_libraries_Fortran(test.f2018_mpi_contracts PRIVATE mpiu90 mod)
    foreach(ranks 1 2 4)
        add_test(NAME f2018_mpi_contracts_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_mpi_contracts> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_mpi_contracts_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
    set(interface_includes)
    foreach(directory IN LISTS MPI_Fortran_INCLUDE_PATH)
        list(APPEND interface_includes "-I${directory}")
    endforeach()
    add_test(NAME f2018_interface_rejections
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_interface_rejections.py
        ${CMAKE_Fortran_COMPILER} "${Fortran_MODULE_DIRECTORY_root}/mod"
        "${Fortran_MODULE_DIRECTORY_root}/mpiu90" ${interface_includes})
    add_executable(test.f2018_system_paths system_paths.f90)
    target_link_libraries_Fortran(test.f2018_system_paths PRIVATE mpiu90 9290 mod)
    add_test(NAME f2018_system_paths
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_system_paths.py
        $<TARGET_FILE:test.f2018_system_paths>)
endif()

if(TARGET mod_csfg)
    add_executable(test.f2018_csfg_reallocation csfg_reallocation.f90)
    target_link_libraries_Fortran(test.f2018_csfg_reallocation PRIVATE mod_csfg mod)
    add_test(NAME f2018_csfg_reallocation COMMAND test.f2018_csfg_reallocation)

    add_executable(test.f2018_nmtejj_csfg nmtejj_probe.f90)
    target_link_libraries_Fortran(test.f2018_nmtejj_csfg PRIVATE rang90_csfg mod)
    if(Python3_Interpreter_FOUND)
        add_test(NAME f2018_nmtejj_csfg_parity
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_nmtejj.py
            $<TARGET_FILE:test.f2018_nmtejj_csfg>)
    endif()
endif()

if(TARGET dvd90mpi)
    add_executable(test.f2018_davidson_mpi davidson_contracts.F90)
    target_compile_definitions(test.f2018_davidson_mpi PRIVATE DAVIDSON_MPI)
    target_link_libraries_Fortran(test.f2018_davidson_mpi PRIVATE dvd90mpi mpiu90 9290 mod)
    target_link_libraries(test.f2018_davidson_mpi PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    foreach(ranks 1 2 4)
        add_test(NAME f2018_davidson_mpi_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_davidson_mpi> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_davidson_mpi_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
endif()
