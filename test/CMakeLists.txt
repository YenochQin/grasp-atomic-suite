add_subdirectory(descriptor_oracle)
# Test executables stay in build/test; application paths use target expressions.
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Let's catch implicit uses of routines in the test prorams and libraries
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=implicit-interface -Werror=implicit-procedure")

# Check lib92 routines -- this currently acts as a simple integration test,
# making sure that everything actually got compiled properly. The test checks
# that the QUAD routine from lib92 actually produces reasonable numbers.
add_executable(test.core_quadrature
    core_quadrature.f90
)
target_link_libraries_Fortran(test.core_quadrature PUBLIC grasp_shared core)
add_test(core_quadrature test.core_quadrature)

add_executable(test.single_subshell_tensor single_subshell_tensor.f90)
target_link_libraries_Fortran(test.single_subshell_tensor PUBLIC spin_angular core grasp_shared)
add_test(NAME single_subshell_tensor COMMAND test.single_subshell_tensor)

add_executable(test.shared_grid_capacity shared_grid_capacity.f90)
target_link_libraries_Fortran(test.shared_grid_capacity PUBLIC grasp_shared)
add_test(NAME shared_grid_capacity COMMAND test.shared_grid_capacity)

if(TARGET grasp_mpi)
    # Calls the production sparse kernels with rank-local column buffers.
    # The test supplies a checked reduction seam, so no MPI job is launched.
    add_executable(test.grasp_mpi_local_sparse grasp_mpi_local_sparse.f90)
    target_link_libraries_Fortran(test.grasp_mpi_local_sparse PUBLIC grasp_mpi grasp_shared core)
    add_test(grasp_mpi_local_sparse test.grasp_mpi_local_sparse)
    set_tests_properties(grasp_mpi_local_sparse PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

add_executable(test.rdensity_natural_orbitals
    ../src/apps/rdensity/natural_orbitals.f90
    rdensity/natural_orbitals_test.f90)
setup_fortran_modules(test.rdensity_natural_orbitals)
target_link_libraries_Fortran(test.rdensity_natural_orbitals PUBLIC grasp_shared core)
add_test(NAME rdensity_natural_orbitals COMMAND test.rdensity_natural_orbitals)

# Prefer the shared workspace environment; standalone installations can select
# their interpreter with -DPython3_EXECUTABLE=...
if(WIN32)
    set(GRASP_GRID_TEST_PYTHON "${PROJECT_SOURCE_DIR}/../graspkit-tools/.venv/Scripts/python.exe")
else()
    set(GRASP_GRID_TEST_PYTHON "${PROJECT_SOURCE_DIR}/../graspkit-tools/.venv/bin/python")
endif()
if(EXISTS "${GRASP_GRID_TEST_PYTHON}" AND NOT Python3_EXECUTABLE)
    set(Python3_EXECUTABLE "${GRASP_GRID_TEST_PYTHON}")
endif()
if(GRASP_STRICT_ACCEPTANCE)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)
else()
    find_package(Python3 COMPONENTS Interpreter)
endif()
if(Python3_Interpreter_FOUND)
    add_test(NAME rlevels_file_output
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rlevels/regression.py"
            $<TARGET_FILE:rlevels> "${PROJECT_SOURCE_DIR}/data")
    set_tests_properties(rlevels_file_output PROPERTIES TIMEOUT 30)
    foreach(variant serial mpi)
        if(variant STREQUAL "mpi" AND NOT TARGET grasp_mpi)
            continue()
        endif()
        set(target test.ris4_cache_${variant})
        add_executable(${target} ris4/cache_io.F90
            "${PROJECT_SOURCE_DIR}/src/apps/isotope/open_isotope_shift_angular_cache.f90"
)
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE core grasp_shared)
        set(cache_mpi_arguments)
        if(variant STREQUAL "mpi")
            target_compile_definitions(${target} PRIVATE CACHE_MPI)
            target_link_libraries_Fortran(${target} PRIVATE grasp_mpi)
            target_link_libraries(${target} PRIVATE MPI::MPI_Fortran)
            set(cache_mpi_arguments "${MPIEXEC_EXECUTABLE}")
        endif()
        add_test(NAME ris4_cache_${variant} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/ris4/check_cache_io.py
            $<TARGET_FILE:${target}> ${cache_mpi_arguments})
        set_tests_properties(ris4_cache_${variant} PROPERTIES TIMEOUT 180)
        if(variant STREQUAL "mpi")
            set_tests_properties(ris4_cache_${variant} PROPERTIES PROCESSORS 2)
        endif()
    endforeach()
    if(TARGET isotope)
        add_test(NAME ris4_cache_reuse_serial
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/ris4/cache_reuse.py"
                $<TARGET_FILE:isotope> "${PROJECT_SOURCE_DIR}/data")
        set_tests_properties(ris4_cache_reuse_serial PROPERTIES TIMEOUT 60)
    endif()
    add_test(NAME fical_mixed_mass_shift_units
        COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/ris4/fical_units.py"
            $<TARGET_FILE:fical>)
    set_tests_properties(fical_mixed_mass_shift_units PROPERTIES TIMEOUT 30 LABELS f2018_properties)

    if(TARGET rdensity)
        add_test(NAME rdensity_serial_regression
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rdensity/regression.py
            --serial $<TARGET_FILE:rdensity> --data ${PROJECT_SOURCE_DIR}/data)
        set_tests_properties(rdensity_serial_regression PROPERTIES TIMEOUT 180)
        if(TARGET rdensity_mpi AND MPIEXEC_EXECUTABLE)
            foreach(ranks 1 2 4)
                add_test(NAME rdensity_mpi_${ranks}
                    COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rdensity/regression.py
                    --serial $<TARGET_FILE:rdensity> --mpi $<TARGET_FILE:rdensity_mpi>
                    --mpiexec ${MPIEXEC_EXECUTABLE} --numproc-flag=${MPIEXEC_NUMPROC_FLAG}
                    --ranks ${ranks} --data ${PROJECT_SOURCE_DIR}/data)
                set_tests_properties(rdensity_mpi_${ranks} PROPERTIES
                    TIMEOUT 180 PROCESSORS ${ranks})
            endforeach()
        endif()
    endif()
    add_test(NAME grasp_grid_patcher
        COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_patch_grasp_grid.py")
    set_tests_properties(grasp_grid_patcher PROPERTIES ENVIRONMENT "FC=${CMAKE_Fortran_COMPILER}")
    add_test(NAME mcp_compare
        COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_mcp_compare.py")
    set_tests_properties(mcp_compare PROPERTIES TIMEOUT 120 SKIP_RETURN_CODE 77
        ENVIRONMENT "FC=${CMAKE_Fortran_COMPILER}")
    add_test(NAME rmcdhf_mpi_sampler
        COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_sample_rmcdhf_mpi.py")
    set_tests_properties(rmcdhf_mpi_sampler PROPERTIES TIMEOUT 30)
    if(TARGET rwfnestimate AND TARGET rangular AND TARGET rci AND TARGET gj)
        add_test(NAME descriptor_v3_calculations_serial COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_v3_calculation_boundaries.py" $<TARGET_FILE_DIR:rangular>)
        if(TARGET rangular_mpi AND TARGET rci_mpi AND TARGET gj_mpi AND MPIEXEC_EXECUTABLE)
            add_test(NAME descriptor_v3_calculations_mpi COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/test_v3_calculation_boundaries.py" $<TARGET_FILE_DIR:rangular>
                --mpiexec "${MPIEXEC_EXECUTABLE}")
            set_tests_properties(descriptor_v3_calculations_mpi PROPERTIES PROCESSORS 2 TIMEOUT 180)
        endif()
    endif()
    if(TARGET rwfnestimate AND TARGET rangular AND TARGET rmcdhf_mem AND TARGET rci)
        add_test(NAME state_workflow_serial COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_state_workflow.py" $<TARGET_FILE_DIR:rwfnestimate>)
        if(TARGET rangular_mpi AND TARGET rmcdhf_mem_mpi AND TARGET rci_mpi AND MPIEXEC_EXECUTABLE)
            add_test(NAME state_workflow_mpi COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/test_state_workflow.py" $<TARGET_FILE_DIR:rwfnestimate>
                --mpiexec "${MPIEXEC_EXECUTABLE}")
            set_tests_properties(state_workflow_mpi PROPERTIES PROCESSORS 2)
        endif()
    endif()
    if(TARGET rwfnestimate AND TARGET rangular_mpi AND TARGET rmcdhf_mpi AND MPIEXEC_EXECUTABLE)
        add_executable(test.indexed_mcp f2018/indexed_mcp_roundtrip.f90
            "${PROJECT_SOURCE_DIR}/src/apps/rangular_mpi/write_mcp_sparsity_summary_mpi.f90")
        setup_fortran_modules(test.indexed_mcp)
        target_link_libraries_Fortran(test.indexed_mcp PRIVATE grasp_mpi core grasp_shared)
        add_test(NAME rmcdhf_mcp_paths COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_rmcdhf_mcp_paths.py" $<TARGET_FILE_DIR:rwfnestimate>
            --mpiexec "${MPIEXEC_EXECUTABLE}" --probe $<TARGET_FILE:test.indexed_mcp>)
        set_tests_properties(rmcdhf_mcp_paths PROPERTIES PROCESSORS 2 TIMEOUT 180)
    endif()
    if(TARGET rwfnestimate AND TARGET rci)
        add_test(NAME rci_restart_serial COMMAND "${Python3_EXECUTABLE}"
            "${CMAKE_CURRENT_SOURCE_DIR}/test_rci_restart.py"
            --rwfnestimate $<TARGET_FILE:rwfnestimate> --rci $<TARGET_FILE:rci>)
        set_tests_properties(rci_restart_serial PROPERTIES TIMEOUT 120 LABELS "f2018_io;f2018_serial")
    endif()
    if(TARGET rwfnestimate AND TARGET rci_mpi AND MPIEXEC_EXECUTABLE)
        foreach(ranks IN ITEMS 1 2 4)
            add_test(NAME rci_restart_mpi_${ranks} COMMAND "${Python3_EXECUTABLE}"
                "${CMAKE_CURRENT_SOURCE_DIR}/test_rci_restart.py"
                --rwfnestimate $<TARGET_FILE:rwfnestimate> --rci $<TARGET_FILE:rci_mpi>
                --mpiexec "${MPIEXEC_EXECUTABLE}" --numproc-flag=${MPIEXEC_NUMPROC_FLAG} --ranks ${ranks})
            set_tests_properties(rci_restart_mpi_${ranks} PROPERTIES TIMEOUT 240 PROCESSORS ${ranks}
                LABELS "f2018_io;f2018_mpi")
        endforeach()
    endif()
    if(TARGET rangular_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME rangular_mpi_descriptors
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rangular/test_mpi_descriptors.py"
            --rangular $<TARGET_FILE:rangular_mpi> --rci $<TARGET_FILE:rci_mpi>
            --mpiexec "${MPIEXEC_EXECUTABLE}"
            --numproc-flag=${MPIEXEC_NUMPROC_FLAG})
        set_tests_properties(rangular_mpi_descriptors PROPERTIES TIMEOUT 180 PROCESSORS 4
            LABELS "f2018_io;f2018_mpi")
    endif()
    if(TARGET rci_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME rci_mpi_breit
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rci_breit.py"
            --rci $<TARGET_FILE:rci_mpi> --mpiexec "${MPIEXEC_EXECUTABLE}"
            --numproc-flag=${MPIEXEC_NUMPROC_FLAG})
        set_tests_properties(rci_mpi_breit PROPERTIES TIMEOUT 120 PROCESSORS 2
            LABELS "f2018_mpi")
    endif()
    if(TARGET csf_descriptor_verify)
        add_test(NAME csf_descriptor_input
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_csf_descriptors.py"
            --rangular $<TARGET_FILE:rangular> --verifier $<TARGET_FILE:csf_descriptor_verify>)
        set_tests_properties(csf_descriptor_input PROPERTIES TIMEOUT 180
            LABELS "f2018_io;f2018_serial"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME rangular_version_comparison
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/test_compare_rangular_versions.py"
            --rangular $<TARGET_FILE:rangular> --verifier $<TARGET_FILE:csf_descriptor_verify>)
        set_tests_properties(rangular_version_comparison PROPERTIES TIMEOUT 180
            LABELS "f2018_io;f2018_serial"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
        add_test(NAME rangular_angular_reuse
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rangular/test_reuse.py"
            --rangular $<TARGET_FILE:rangular> --verifier $<TARGET_FILE:csf_descriptor_verify>)
        set_tests_properties(rangular_angular_reuse PROPERTIES TIMEOUT 180
            LABELS "f2018_io;f2018_serial"
            ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
    endif()
endif()

add_subdirectory(jj2lsj)
add_subdirectory(standard_grasp)
add_subdirectory(rwfnestimate)
add_subdirectory(f2018)

# Analytical EAL weights exercise the former absent-JCUPA dependency.
foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(TARGET ${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
        add_executable(test.descriptor_weights_${variant} descriptor_weights.f90
            "${app}/read_al_level_weights.f90"
            "${app}/basis_csf_two_j_plus_one.f90")
        setup_fortran_modules(test.descriptor_weights_${variant})
        target_link_libraries_Fortran(test.descriptor_weights_${variant} PRIVATE core grasp_shared)
        add_test(NAME descriptor_weights_${variant} COMMAND test.descriptor_weights_${variant})
    endif()
endforeach()
if(TARGET csf_descriptor_verify AND Python3_Interpreter_FOUND)
    add_test(NAME native_descriptor_integration COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_native_descriptors.py --bin $<TARGET_FILE_DIR:rangular>)
endif()

foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(TARGET ${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
        set(target test.descriptor_input_${variant})
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            set(loader load_csf_basis_data_mpi)
            set(angular rangular_mpi)
        else()
            set(loader load_csf_basis_data)
            set(angular rangular)
        endif()
        add_executable(${target} descriptor_rmcdhf.F90 "${app}/${loader}.f90"
            "${app}/connect_mcp_files.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE core grasp_shared)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE DESCRIPTOR_MPI)
            target_link_libraries_Fortran(${target} PRIVATE grasp_mpi)
            set(mpi_arguments --mpiexec "${MPIEXEC_EXECUTABLE}")
        endif()
        if(NOT variant MATCHES "_mpi$")
            target_sources(${target} PRIVATE "${app}/serial_mpi_context.f90")
        endif()
        if(variant MATCHES "mem")
            target_sources(${target} PRIVATE "${app}/cache_mcp_files.f90" "${app}/rmcdhf_mcp_memory.f90")
            if(variant STREQUAL "rmcdhf_mem")
                target_sources(${target} PRIVATE "${app}/query_total_system_memory.f90"
                    "${app}/query_available_system_memory.f90")
            endif()
        endif()
        add_test(NAME descriptor_input_${variant} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_rmcdhf_descriptors.py
            --angular $<TARGET_FILE:${angular}> --loader $<TARGET_FILE:${target}> ${mpi_arguments})
        set_tests_properties(descriptor_input_${variant} PROPERTIES TIMEOUT 180)
    endif()
endforeach()

add_executable(test.descriptor_lifecycle descriptor_lifecycle.f90)
target_link_libraries_Fortran(test.descriptor_lifecycle PRIVATE core grasp_shared)

# Also check the production merge source under the F2018 line limit in Release builds.
add_test(NAME descriptor_merge_f2018 COMMAND "${CMAKE_Fortran_COMPILER}"
    -std=f2018 -ffree-line-length-132 -Werror=line-truncation -fsyntax-only
    -I$<TARGET_PROPERTY:grasp_shared,Fortran_MODULE_DIRECTORY>
    "${PROJECT_SOURCE_DIR}/src/lib/core/csf_basis_merging.f90")

if(Python3_Interpreter_FOUND)
    if(TARGET gj AND TARGET rci)
        add_test(NAME lande_factors_serial COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_lande_factors.py --gj $<TARGET_FILE:gj> --rci $<TARGET_FILE:rci>)
        set_tests_properties(lande_factors_serial PROPERTIES TIMEOUT 180)
    endif()
    if(TARGET gj_mpi AND TARGET rci_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME lande_factors_mpi COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_lande_factors.py --gj $<TARGET_FILE:gj_mpi> --rci $<TARGET_FILE:rci_mpi>
            --mpiexec "${MPIEXEC_EXECUTABLE}")
        set_tests_properties(lande_factors_mpi PROPERTIES TIMEOUT 180 PROCESSORS 2)
    endif()
    if(TARGET isotope AND TARGET rci AND TARGET fical)
        add_test(NAME isotope_shifts_serial COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_isotope_shifts.py --ris $<TARGET_FILE:isotope> --rci $<TARGET_FILE:rci>
            --fical $<TARGET_FILE:fical>)
        set_tests_properties(isotope_shifts_serial PROPERTIES TIMEOUT 180)
    endif()
    if(TARGET isotope_mpi AND TARGET rci_mpi AND TARGET fical AND MPIEXEC_EXECUTABLE)
        add_test(NAME isotope_shifts_mpi COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_isotope_shifts.py --ris $<TARGET_FILE:isotope_mpi> --rci $<TARGET_FILE:rci_mpi>
            --fical $<TARGET_FILE:fical> --mpiexec "${MPIEXEC_EXECUTABLE}")
        set_tests_properties(isotope_shifts_mpi PROPERTIES TIMEOUT 180 PROCESSORS 2)
    endif()
    if(TARGET hyperfine AND TARGET rci)
        add_test(NAME hyperfine_constants_serial COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_hyperfine_constants.py --rhfs $<TARGET_FILE:hyperfine> --rci $<TARGET_FILE:rci>)
        set_tests_properties(hyperfine_constants_serial PROPERTIES TIMEOUT 180)
    endif()
    if(TARGET hyperfine_mpi AND TARGET rci_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME hyperfine_constants_mpi COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_hyperfine_constants.py --rhfs $<TARGET_FILE:hyperfine_mpi> --rci $<TARGET_FILE:rci_mpi>
            --mpiexec "${MPIEXEC_EXECUTABLE}")
        set_tests_properties(hyperfine_constants_mpi PROPERTIES TIMEOUT 180 PROCESSORS 2)
    endif()
    if(TARGET biorthtransform)
        add_test(NAME biorthogonal_states_serial COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_biorthogonal_states.py $<TARGET_FILE:biorthtransform>)
        set_tests_properties(biorthogonal_states_serial PROPERTIES TIMEOUT 180)
    endif()
    if(TARGET biorthtransform_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME biorthogonal_states_mpi COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_biorthogonal_states.py $<TARGET_FILE:biorthtransform_mpi>
            --mpiexec "${MPIEXEC_EXECUTABLE}")
        set_tests_properties(biorthogonal_states_mpi PROPERTIES TIMEOUT 180 PROCESSORS 2)
    endif()
    if(TARGET biorthtransform AND TARGET rtransition AND TARGET rtransition_phase)
        add_test(NAME radiative_transitions_serial COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_radiative_transitions.py $<TARGET_FILE:rtransition>)
        add_test(NAME radiative_transitions_phase COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_radiative_transitions.py $<TARGET_FILE:rtransition_phase>)
        set_tests_properties(radiative_transitions_serial radiative_transitions_phase PROPERTIES TIMEOUT 180)
    endif()
    if(TARGET biorthtransform_mpi AND TARGET rtransition_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME radiative_transitions_mpi COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_radiative_transitions.py $<TARGET_FILE:rtransition_mpi>
            --mpiexec "${MPIEXEC_EXECUTABLE}")
        set_tests_properties(radiative_transitions_mpi PROPERTIES TIMEOUT 180 PROCESSORS 2)
    endif()
endif()

if(TARGET biorthtransform AND TARGET rtransition AND TARGET rtransition_phase AND Python3_Interpreter_FOUND)
    add_test(NAME descriptor_dual_basis_pipeline COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_descriptor_pipeline.py --bin $<TARGET_FILE_DIR:biorthtransform>)
    set_tests_properties(descriptor_dual_basis_pipeline PROPERTIES TIMEOUT 180)
    if(TARGET biorthtransform_mpi AND TARGET rtransition_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME descriptor_dual_basis_pipeline_mpi COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_descriptor_pipeline.py --bin $<TARGET_FILE_DIR:biorthtransform_mpi>
            --mpiexec "${MPIEXEC_EXECUTABLE}" --ci)
        set_tests_properties(descriptor_dual_basis_pipeline_mpi PROPERTIES TIMEOUT 180 PROCESSORS 2)
    endif()
endif()
if(TARGET grasp_mpi)
    add_executable(test.descriptor_collectives descriptor_collectives.f90)
    target_link_libraries_Fortran(test.descriptor_collectives PRIVATE grasp_mpi core grasp_shared)
    foreach(ranks 1 2 4)
        add_test(NAME descriptor_collectives_mpi_${ranks}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
            $<TARGET_FILE:test.descriptor_collectives> valid)
        set_tests_properties(descriptor_collectives_mpi_${ranks} PROPERTIES TIMEOUT 30 PROCESSORS ${ranks})
    endforeach()
    foreach(scenario error shape)
        add_test(NAME descriptor_collectives_mpi_${scenario}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2
            $<TARGET_FILE:test.descriptor_collectives> ${scenario})
        set_tests_properties(descriptor_collectives_mpi_${scenario} PROPERTIES TIMEOUT 30 PROCESSORS 2 WILL_FAIL TRUE)
    endforeach()
endif()

foreach(variant rtransition rtransition_phase rtransition_mpi)
    if(TARGET ${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/apps/${variant}")
        set(target test.descriptor_merge_${variant})
        add_executable(${target} descriptor_merge.F90
            "${app}/merge_transition_csf_bases.f90"
            "${app}/load_initial_transition_basis.f90"
            "${app}/load_final_transition_basis.f90"
            "${app}/merge_transition_basis_records.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE core grasp_shared)
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE DESCRIPTOR_MPI)
            target_link_libraries_Fortran(${target} PRIVATE grasp_mpi)
            set(mpi_arguments --mpiexec "${MPIEXEC_EXECUTABLE}")
        endif()
        add_test(NAME descriptor_merge_${variant} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_descriptor_merge.py $<TARGET_FILE:${target}> ${mpi_arguments})
        set_tests_properties(descriptor_merge_${variant} PROPERTIES TIMEOUT 60)
    endif()
endforeach()

set(RCSFS_TEST_CLI "$ENV{RCSFS_TEST_CLI}" CACHE FILEPATH "V3-capable rcsfs CLI used only to prepare test fixtures")
function(set_native_fixture_environment directory)
    get_property(tests DIRECTORY "${directory}" PROPERTY TESTS)
    if(tests AND RCSFS_TEST_CLI)
        set_property(TEST ${tests} DIRECTORY "${directory}" APPEND PROPERTY ENVIRONMENT "RCSFS_TEST_CLI=${RCSFS_TEST_CLI}")
    endif()
    get_property(children DIRECTORY "${directory}" PROPERTY SUBDIRECTORIES)
    foreach(child IN LISTS children)
        set_native_fixture_environment("${child}")
    endforeach()
endfunction()
set_native_fixture_environment("${CMAKE_CURRENT_SOURCE_DIR}")
