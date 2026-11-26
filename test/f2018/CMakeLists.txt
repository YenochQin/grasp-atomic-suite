add_executable(test.f2018_library_modules library_modules.f90)
target_link_libraries_Fortran(test.f2018_library_modules PRIVATE core grasp_shared)
add_test(NAME f2018_library_modules COMMAND test.f2018_library_modules)

add_executable(test.f2018_breit_procedure_hooks breit_procedure_hooks.f90)
target_link_libraries_Fortran(test.f2018_breit_procedure_hooks PRIVATE breit_angular_support core grasp_shared)
add_test(NAME f2018_breit_procedure_hooks COMMAND test.f2018_breit_procedure_hooks)

add_executable(test.f2018_davidson davidson_contracts.F90)
target_link_libraries_Fortran(test.f2018_davidson PRIVATE davidson_eigensolver core grasp_shared)
target_link_libraries(test.f2018_davidson PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME f2018_davidson_contracts COMMAND test.f2018_davidson)
add_test(NAME f2018_davidson_early COMMAND test.f2018_davidson early)
set_tests_properties(f2018_davidson_contracts f2018_davidson_early PROPERTIES TIMEOUT 30)
foreach(variant rci rci_mpi)
    set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
    set(target test.f2018_rk_count_${variant})
    add_executable(${target} rk_count.F90 rk_count_doubles.f90
        "${app}/tabulate_slater_integrals.f90" "${app}/satisfies_coulomb_triangle.f90")
    set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
    target_link_libraries_Fortran(${target} PRIVATE core grasp_shared)
    add_test(NAME f2018_rk_count_${variant} COMMAND ${target})
    set_tests_properties(f2018_rk_count_${variant} PROPERTIES TIMEOUT 30)
endforeach()
add_executable(test.f2018_nmtejj nmtejj_probe.f90)
add_executable(test.f2018_linear_algebra linear_algebra_contracts.f90)
target_link_libraries_Fortran(test.f2018_linear_algebra PRIVATE grasp_shared)
target_link_libraries(test.f2018_linear_algebra PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME f2018_linear_algebra_contracts COMMAND test.f2018_linear_algebra)
foreach(variant rbiotransform rbiotransform_mpi)
    add_executable(test.bndinv_${variant} bndinv_precision.f90
        "${PROJECT_SOURCE_DIR}/src/apps/${variant}/invert_banded_matrix.f90")
    setup_fortran_modules(test.bndinv_${variant})
    add_test(NAME f2018_bndinv_${variant} COMMAND test.bndinv_${variant})
endforeach()
if(Python3_Interpreter_FOUND)
    add_test(NAME f2018_library_source_style COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_fortran_source_style.py lib)
    add_test(NAME f2018_application_source_style COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_fortran_source_style.py apps)
    add_test(NAME f2018_numeric_literals COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_numeric_literals.py ${CMAKE_Fortran_COMPILER})
    add_test(NAME f2018_bndinv_interfaces COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_bndinv_interface.py ${CMAKE_Fortran_COMPILER})
    add_test(NAME f2018_numeric_provenance COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/test/npol/check_upstream_tamper.py)
    add_test(NAME f2018_qed_boundary COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_qed_boundary.py ${CMAKE_Fortran_COMPILER}
        "${Fortran_MODULE_DIRECTORY_root}/grasp_shared" $<TARGET_FILE:grasp_shared>)
    set_tests_properties(f2018_numeric_literals f2018_qed_boundary PROPERTIES TIMEOUT 60)
endif()
add_executable(test.f2018_storage_bytes storage_bytes.f90)
target_link_libraries_Fortran(test.f2018_storage_bytes PRIVATE grasp_shared)
add_test(NAME f2018_storage_bytes COMMAND test.f2018_storage_bytes)
target_link_libraries_Fortran(test.f2018_nmtejj PRIVATE spin_angular grasp_shared)
if(Python3_Interpreter_FOUND)
    add_test(NAME f2018_nmtejj_parity
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_nmtejj.py
        $<TARGET_FILE:test.f2018_nmtejj>)
endif()

if(Python3_Interpreter_FOUND)
    add_executable(test.f2018_interp interp_contracts.f90
        "${PROJECT_SOURCE_DIR}/src/lib/core/aitken_interpolation.f90")
    target_link_libraries_Fortran(test.f2018_interp PRIVATE grasp_shared)
    target_compile_options(test.f2018_interp PRIVATE -finit-integer=999999 -fcheck=all)
    add_test(NAME f2018_interp_contracts COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_interp.py $<TARGET_FILE:test.f2018_interp>)
    set_tests_properties(f2018_interp_contracts PROPERTIES TIMEOUT 60)
    foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
        set(target test.f2018_weights_${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
        add_executable(${target} level_weights.f90 level_weights_doubles.f90 "${app}/read_al_level_weights.f90")
        set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
        target_link_libraries_Fortran(${target} PRIVATE core grasp_shared)
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
foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(variant MATCHES "_mpi$" AND NOT TARGET grasp_mpi)
        continue()
    endif()
    set(directory ${variant})
    set(app "${PROJECT_SOURCE_DIR}/src/apps/${directory}")
    set(matrix solve_hamiltonian_blocks)
    if(variant MATCHES "_mpi$")
        set(matrix solve_hamiltonian_blocks_mpi)
    endif()
    set(target test.f2018_ci_${variant})
    add_executable(${target} ci_damping.F90 ci_damping_doubles.F90 "${app}/${matrix}.f90")
    set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
    target_link_libraries_Fortran(${target} PRIVATE core grasp_shared)
    target_link_libraries(${target} PRIVATE ${BLAS_LIBRARIES})
    target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace)
    if(variant MATCHES "mem")
        target_compile_definitions(${target} PRIVATE CI_MEM)
        target_sources(${target} PRIVATE "${app}/rmcdhf_mcp_memory.f90")
    endif()
    if(variant STREQUAL "rmcdhf_mpi")
        target_compile_definitions(${target} PRIVATE CI_INDEXED)
        target_sources(${target} PRIVATE "${app}/rmcdhf_mcp_io.f90" "${app}/mcp_density_contractions.f90"
            "${app}/csf_density_matrix_element.f90")
    endif()
    if(variant MATCHES "_mpi$")
        target_compile_definitions(${target} PRIVATE CI_MPI)
        target_link_libraries_Fortran(${target} PRIVATE grasp_mpi)
        foreach(ranks 1 2 4)
            add_test(NAME f2018_ci_damping_${variant}_${ranks}
                COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${target}> ${MPIEXEC_POSTFLAGS})
            set(fixture_dir "${CMAKE_CURRENT_BINARY_DIR}/ci-damping-${variant}-${ranks}")
            file(MAKE_DIRECTORY "${fixture_dir}")
            set_tests_properties(f2018_ci_damping_${variant}_${ranks} PROPERTIES
                TIMEOUT 60 PROCESSORS ${ranks} WORKING_DIRECTORY "${fixture_dir}")
        endforeach()
    else()
        target_sources(${target} PRIVATE "${app}/serial_mpi_context.f90")
        add_test(NAME f2018_ci_damping_${variant} COMMAND ${target})
        set(fixture_dir "${CMAKE_CURRENT_BINARY_DIR}/ci-damping-${variant}")
        file(MAKE_DIRECTORY "${fixture_dir}")
        set_tests_properties(f2018_ci_damping_${variant} PROPERTIES TIMEOUT 30 WORKING_DIRECTORY "${fixture_dir}")
    endif()
endforeach()

if(TARGET grasp_mpi AND Python3_Interpreter_FOUND)
    add_test(NAME f2018_worker_output COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_worker_output.py ${CMAKE_Fortran_COMPILER}
        ${MPI_Fortran_COMPILER} ${MPIEXEC_EXECUTABLE})
    set_tests_properties(f2018_worker_output PROPERTIES TIMEOUT 90 PROCESSORS 2)
    set(angular_app "${PROJECT_SOURCE_DIR}/src/apps/rangular_mpi")
    add_executable(test.f2018_sparse_empty_rank sparse_empty_rank.f90
        "${angular_app}/write_mcp_sparsity.f90" "${angular_app}/write_mcp_sparsity_summary_mpi.f90")
    set_property(TARGET test.f2018_sparse_empty_rank PROPERTY Fortran_MODULE_DIRECTORY
        "${Fortran_MODULE_DIRECTORY_root}/test.f2018_sparse_empty_rank")
    target_link_libraries_Fortran(test.f2018_sparse_empty_rank PRIVATE core grasp_shared grasp_mpi)
    # Poison uninitialized local integers only in this diagnostic fixture.
    # A zero default would hide SETSDA's former empty-rank IC access.
    target_compile_options(test.f2018_sparse_empty_rank PRIVATE -finit-integer=-999999 -fcheck=all)
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sparse-empty-rank")
    add_test(NAME f2018_sparse_empty_rank COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_sparse_empty_rank> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(f2018_sparse_empty_rank PROPERTIES TIMEOUT 30 PROCESSORS 2
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sparse-empty-rank")
    file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sparse-empty-rank-4")
    add_test(NAME f2018_sparse_empty_rank_4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
        ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_sparse_empty_rank> ${MPIEXEC_POSTFLAGS})
    set_tests_properties(f2018_sparse_empty_rank_4 PROPERTIES TIMEOUT 30 PROCESSORS 4
        WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/sparse-empty-rank-4")
    add_executable(test.f2018_mpi_contracts mpi_contracts.f90)
    target_link_libraries_Fortran(test.f2018_mpi_contracts PRIVATE grasp_mpi grasp_shared)
    foreach(ranks 1 2 4)
        add_test(NAME f2018_mpi_contracts_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:test.f2018_mpi_contracts> ${MPIEXEC_POSTFLAGS})
        set_tests_properties(f2018_mpi_contracts_${ranks} PROPERTIES TIMEOUT 60 PROCESSORS ${ranks})
    endforeach()
    foreach(ranks 2 4)
        add_test(NAME f2018_mpi_count_overflow_${ranks}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mpi_count_failure.py
            ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks} ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:test.f2018_mpi_contracts> ${MPIEXEC_POSTFLAGS} overflow)
        set_tests_properties(f2018_mpi_count_overflow_${ranks} PROPERTIES TIMEOUT 30 PROCESSORS ${ranks})
    endforeach()
    set(interface_includes)
    foreach(directory IN LISTS MPI_Fortran_INCLUDE_PATH)
        list(APPEND interface_includes "-I${directory}")
    endforeach()
    add_test(NAME f2018_interface_rejections
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_interface_rejections.py
        ${CMAKE_Fortran_COMPILER} "${Fortran_MODULE_DIRECTORY_root}/grasp_shared"
        "${Fortran_MODULE_DIRECTORY_root}/grasp_mpi" ${interface_includes})
    add_executable(test.f2018_system_paths system_paths.f90)
    target_link_libraries_Fortran(test.f2018_system_paths PRIVATE grasp_mpi core grasp_shared)
    add_test(NAME f2018_system_paths
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_system_paths.py
        $<TARGET_FILE:test.f2018_system_paths>)
endif()

# Input failure tests call the production routines with valid, empty and
# truncated records; a single damaged MPI rank must abort without a barrier.
if(Python3_Interpreter_FOUND)
    add_test(NAME f2018_input_precision
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_input_precision.py
        ${CMAKE_Fortran_COMPILER} "${Fortran_MODULE_DIRECTORY_root}/grasp_shared")
    set_tests_properties(f2018_input_precision PROPERTIES TIMEOUT 60)
    foreach(variant rmcdhf rmcdhf_mem)
        if(NOT TARGET ${variant})
            continue()
        endif()
        set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
        set(target test.mcp_io_${variant})
        add_executable(${target} mcp_reader_io.F90 "${app}/connect_mcp_files.f90")
        set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
        target_link_libraries_Fortran(${target} PRIVATE grasp_shared core)
        target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace -finit-integer=0)
        target_sources(${target} PRIVATE "${app}/serial_mpi_context.f90")
        if(variant MATCHES "mem")
            target_compile_definitions(${target} PRIVATE MCP_MEM)
            target_sources(${target} PRIVATE mcp_memory_queries.f90
                "${app}/cache_mcp_files.f90" "${app}/rmcdhf_mcp_memory.f90")
        else()
            target_sources(${target} PRIVATE mcp_angular_doubles.f90 "${app}/assemble_orbital_potential_coefficients.f90")
        endif()
        add_test(NAME f2018_mcp_io_${variant}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mcp_io.py
            $<TARGET_FILE:${target}> --reader)
        set_tests_properties(f2018_mcp_io_${variant} PROPERTIES TIMEOUT 240)
        if(variant MATCHES "mem")
            set(target test.mcp_memory_${variant})
            add_executable(${target} mcp_memory_io.F90 mcp_memory_io_injection.f90 mcp_memory_queries.f90
                "${app}/cache_mcp_files.f90" "${app}/rmcdhf_mcp_memory.f90")
            set_property(TARGET ${target} PROPERTY Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/${target}")
            target_link_libraries_Fortran(${target} PRIVATE grasp_shared)
            target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace -finit-integer=0)
            target_sources(${target} PRIVATE "${app}/serial_mpi_context.f90")
            add_test(NAME f2018_mcp_memory_${variant}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mcp_io.py
                $<TARGET_FILE:${target}>)
            set_tests_properties(f2018_mcp_memory_${variant} PROPERTIES TIMEOUT 180)
        endif()
    endforeach()
    foreach(variant rbiotransform rbiotransform_mpi)
        if(NOT TARGET ${variant})
            continue()
        endif()
        set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
        set(target test.csf_record_${variant})
        add_executable(${target} csf_record_io.F90
            "${PROJECT_SOURCE_DIR}/test/descriptor_oracle/${variant}/lodcslBio.f90"
            "${PROJECT_SOURCE_DIR}/test/descriptor_oracle/${variant}/lodcslBio_I.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE csf_text_oracle grasp_shared core)
        target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace)
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE MCP_MPI)
            target_link_libraries_Fortran(${target} PRIVATE grasp_mpi)
            set(mpi_arguments --mpiexec "${MPIEXEC_EXECUTABLE}")
        endif()
        add_test(NAME f2018_csf_record_${variant}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_csf_record_io.py
            $<TARGET_FILE:${target}> ${mpi_arguments})
        set_tests_properties(f2018_csf_record_${variant} PROPERTIES TIMEOUT 90)
        if(variant MATCHES "_mpi$")
            set_tests_properties(f2018_csf_record_${variant} PROPERTIES PROCESSORS 2)
        endif()
    endforeach()
    if(TARGET grasp_mpi)
        add_executable(test.cpath_io cpath_io.f90)
        target_link_libraries_Fortran(test.cpath_io PRIVATE grasp_mpi grasp_shared)
        add_test(NAME f2018_cpath_io_mpi
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_cpath_io.py
            "${MPIEXEC_EXECUTABLE}" $<TARGET_FILE:test.cpath_io>)
        set_tests_properties(f2018_cpath_io_mpi PROPERTIES TIMEOUT 60 PROCESSORS 4)
    endif()
endif()
