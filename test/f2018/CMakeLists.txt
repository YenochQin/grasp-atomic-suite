add_executable(test.f2018_davidson davidson_contracts.F90)
target_link_libraries_Fortran(test.f2018_davidson PRIVATE dvd90 9290 mod)
target_link_libraries(test.f2018_davidson PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME f2018_davidson_contracts COMMAND test.f2018_davidson)
add_test(NAME f2018_davidson_early COMMAND test.f2018_davidson early)
set_tests_properties(f2018_davidson_contracts f2018_davidson_early PROPERTIES TIMEOUT 30)
foreach(variant rci90 rci90_mpi)
    set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
    set(target test.f2018_rk_count_${variant})
    add_executable(${target} rk_count.F90 "${app}/genintrk.f90" "${app}/genintrk_I.f90"
        "${app}/triangrk.f90" "${app}/triangrk_I.f90")
    target_link_libraries_Fortran(${target} PRIVATE 9290 mod)
    add_test(NAME f2018_rk_count_${variant} COMMAND ${target})
    set_tests_properties(f2018_rk_count_${variant} PROPERTIES TIMEOUT 30)
endforeach()
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
    add_executable(test.f2018_interp interp_contracts.f90
        "${PROJECT_SOURCE_DIR}/src/lib/lib9290/interp.f90"
        "${PROJECT_SOURCE_DIR}/src/lib/lib9290/interp_I.f90")
    target_link_libraries_Fortran(test.f2018_interp PRIVATE mod)
    target_compile_options(test.f2018_interp PRIVATE -finit-integer=999999 -fcheck=all)
    add_test(NAME f2018_interp_contracts COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_interp.py $<TARGET_FILE:test.f2018_interp>)
    set_tests_properties(f2018_interp_contracts PROPERTIES TIMEOUT 60)
    foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi rmcdhf90_mem_strict_mpi rmcdhf90_mem_diag_mpi rmcdhf90_csfg_mpi)
        set(target test.f2018_weights_${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
        add_executable(${target} level_weights.f90 "${app}/getaldwt.f90" "${app}/getaldwt_I.f90")
        target_link_libraries_Fortran(${target} PRIVATE 9290 mod)
        # Diagnostic poison exposes the former uninitialized read extent;
        # no production target receives default-initialization options.
        target_compile_options(${target} PRIVATE -finit-integer=-999999 -fcheck=all)
        add_test(NAME f2018_weights_${variant} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_level_weights.py $<TARGET_FILE:${target}>)
        set_tests_properties(f2018_weights_${variant} PROPERTIES TIMEOUT 60)
    endforeach()
    foreach(driver acceptance first_sweep matrix_warnings)
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
    set(strict_app "${PROJECT_SOURCE_DIR}/src/appl/rmcdhf90_mem_strict_mpi")
    add_executable(test.f2018_strict_orbital_consistency strict_orbital_consistency.f90)
    foreach(source improvmpi consis dampor)
        target_sources(test.f2018_strict_orbital_consistency PRIVATE "${strict_app}/${source}.f90")
    endforeach()
    foreach(interface improvmpi cofpotmpi defcor solve consis dampck dampor orthy setlagmpi)
        target_sources(test.f2018_strict_orbital_consistency PRIVATE "${strict_app}/${interface}_I.f90")
    endforeach()
    target_link_libraries_Fortran(test.f2018_strict_orbital_consistency PRIVATE mpiu90 9290 mod)
    add_test(NAME f2018_strict_orbital_consistency COMMAND test.f2018_strict_orbital_consistency)
    set_tests_properties(f2018_strict_orbital_consistency PROPERTIES TIMEOUT 30)
    set(angular_app "${PROJECT_SOURCE_DIR}/src/appl/rangular90_mpi")
    add_executable(test.f2018_sparse_empty_rank sparse_empty_rank.f90
        "${angular_app}/setsda.f90" "${angular_app}/setsda_I.f90" "${angular_app}/outsdampi_I.f90")
    target_link_libraries_Fortran(test.f2018_sparse_empty_rank PRIVATE 9290 mod)
    # Poison uninitialized local integers only in this diagnostic fixture.
    # A zero default would hide SETSDA's former empty-rank IC access.
    target_compile_options(test.f2018_sparse_empty_rank PRIVATE -finit-integer=-999999 -fcheck=all)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sparse-empty-rank")
    add_test(NAME f2018_sparse_empty_rank COMMAND test.f2018_sparse_empty_rank)
    set_tests_properties(f2018_sparse_empty_rank PROPERTIES TIMEOUT 30
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sparse-empty-rank")
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
    set(app "${PROJECT_SOURCE_DIR}/src/lib/libme_csfg")
    add_executable(test.f2018_rk_count_csfg rk_count.F90 "${app}/genintrk.f90" "${app}/genintrk_I.f90"
        "${app}/triangrk.f90" "${app}/triangrk_I.f90")
    target_compile_definitions(test.f2018_rk_count_csfg PRIVATE RK_CSFG)
    target_link_libraries_Fortran(test.f2018_rk_count_csfg PRIVATE mod_csfg 9290 mod)
    add_test(NAME f2018_rk_count_csfg COMMAND test.f2018_rk_count_csfg)
    set_tests_properties(f2018_rk_count_csfg PROPERTIES TIMEOUT 30)
    set(app "${PROJECT_SOURCE_DIR}/src/appl/rmcdhf90_csfg_mpi")
    add_executable(test.f2018_csfg_dacon csfg_dacon.f90 "${app}/dacon.f90" "${app}/dacon_I.f90")
    target_link_libraries_Fortran(test.f2018_csfg_dacon PRIVATE mod_csfg 9290 mod)
    add_test(NAME f2018_csfg_dacon COMMAND test.f2018_csfg_dacon)
    set_tests_properties(f2018_csfg_dacon PROPERTIES TIMEOUT 30)
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
    add_executable(test.f2018_dgemvmpi dgemvmpi_contracts.f90)
    target_link_libraries_Fortran(test.f2018_dgemvmpi PRIVATE dvd90mpi mpiu90 9290 mod)
    target_link_libraries(test.f2018_dgemvmpi PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    foreach(ranks 1 2 4)
        add_test(NAME f2018_dgemvmpi_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_dgemvmpi> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_dgemvmpi_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
    set(app "${PROJECT_SOURCE_DIR}/src/appl/rmcdhf90_csfg_mpi")
    add_executable(test.f2018_csfg_debug_mpi csfg_debug_mpi.f90 "${app}/setdbgmpi.f90"
        "${app}/setdbgmpi_I.f90" "${app}/setdbg_I.f90")
    target_link_libraries_Fortran(test.f2018_csfg_debug_mpi PRIVATE mod_csfg mpiu90 9290 mod)
    foreach(ranks 1 2 4)
        add_test(NAME f2018_csfg_debug_mpi_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_csfg_debug_mpi> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_csfg_debug_mpi_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
    add_executable(test.f2018_davidson_mpi davidson_contracts.F90)
    target_compile_definitions(test.f2018_davidson_mpi PRIVATE DAVIDSON_MPI)
    target_link_libraries_Fortran(test.f2018_davidson_mpi PRIVATE dvd90mpi mpiu90 9290 mod)
    target_link_libraries(test.f2018_davidson_mpi PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    foreach(ranks 1 2 4)
        add_test(NAME f2018_davidson_mpi_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_davidson_mpi> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_davidson_mpi_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
        add_test(NAME f2018_davidson_early_mpi_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_davidson_mpi> early ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_davidson_early_mpi_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
endif()
