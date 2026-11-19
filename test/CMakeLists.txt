add_subdirectory(descriptor_oracle)
# Test executables stay in build/test; application paths use target expressions.
unset(CMAKE_RUNTIME_OUTPUT_DIRECTORY)

# Let's catch implicit uses of routines in the test prorams and libraries
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Werror=implicit-interface -Werror=implicit-procedure")

# Check lib92 routines -- this currently acts as a simple integration test,
# making sure that everything actually got compiled properly. The test checks
# that the QUAD routine from lib92 actually produces reasonable numbers.
add_executable(test.lib9290_quad
    lib9290_quad.f90
)
target_link_libraries_Fortran(test.lib9290_quad PUBLIC mod 9290)
add_test(lib9290_quad test.lib9290_quad)

add_executable(test.shared_grid_capacity shared_grid_capacity.f90)
target_link_libraries_Fortran(test.shared_grid_capacity PUBLIC mod)
add_test(NAME shared_grid_capacity COMMAND test.shared_grid_capacity)

if(TARGET mpiu90)
    # Calls the production sparse kernels with rank-local column buffers.
    # The test supplies a checked reduction seam, so no MPI job is launched.
    add_executable(test.mpi90_local_sparse mpi90_local_sparse.f90)
    target_link_libraries_Fortran(test.mpi90_local_sparse PUBLIC mpiu90 mod 9290)
    add_test(mpi90_local_sparse test.mpi90_local_sparse)
    set_tests_properties(mpi90_local_sparse PROPERTIES
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
endif()

add_executable(test.rdensity_natural_orbitals
    ../src/appl/rdensity/natural_orbitals.f90
    rdensity/natural_orbitals_test.f90)
setup_fortran_modules(test.rdensity_natural_orbitals)
target_link_libraries_Fortran(test.rdensity_natural_orbitals PUBLIC mod 9290)
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
        if(variant STREQUAL "mpi" AND NOT TARGET mpiu90)
            continue()
        endif()
        set(target test.ris4_cache_${variant})
        add_executable(${target} ris4/cache_io.F90
            "${PROJECT_SOURCE_DIR}/src/appl/ris4/angdata.f90"
            "${PROJECT_SOURCE_DIR}/src/appl/ris4/angdata_I.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE mod)
        set(cache_mpi_arguments)
        if(variant STREQUAL "mpi")
            target_compile_definitions(${target} PRIVATE CACHE_MPI)
            target_link_libraries_Fortran(${target} PRIVATE mpiu90)
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
    set_tests_properties(mcp_compare PROPERTIES TIMEOUT 120 SKIP_RETURN_CODE 77)
    add_test(NAME rmcdhf_mpi_sampler
        COMMAND "${Python3_EXECUTABLE}"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_sample_rmcdhf_mpi.py")
    set_tests_properties(rmcdhf_mpi_sampler PROPERTIES TIMEOUT 30)
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

# Independent synthetic NP tests; no GRASP preparation tools or atomic data.
add_subdirectory(npol)
add_subdirectory(jj2lsj)
add_subdirectory(standard_grasp)
add_subdirectory(rwfnestimate)
add_subdirectory(f2018)

# Analytical EAL weights exercise the former absent-JCUPA dependency.
foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(TARGET ${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
        add_executable(test.descriptor_weights_${variant} descriptor_weights.f90
            "${app}/getaldwt.f90" "${app}/getaldwt_I.f90"
            "${app}/itjpo.f90" "${app}/itjpo_I.f90")
        setup_fortran_modules(test.descriptor_weights_${variant})
        target_link_libraries_Fortran(test.descriptor_weights_${variant} PRIVATE 9290 mod)
        add_test(NAME descriptor_weights_${variant} COMMAND test.descriptor_weights_${variant})
    endif()
endforeach()
if(TARGET csf_descriptor_verify AND Python3_Interpreter_FOUND)
    add_test(NAME native_descriptor_integration COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_native_descriptors.py --bin $<TARGET_FILE_DIR:rangular>)
endif()

foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(TARGET ${variant})
        set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
        set(target test.descriptor_input_${variant})
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            set(loader setcslmpi)
            set(angular rangular_mpi)
        else()
            set(loader setcsl)
            set(angular rangular)
        endif()
        add_executable(${target} descriptor_rmcdhf.F90 "${app}/${loader}.f90"
            "${app}/${loader}_I.f90" "${app}/setmcp.f90" "${app}/setmcp_I.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE 9290 mod)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE DESCRIPTOR_MPI)
            target_link_libraries_Fortran(${target} PRIVATE mpiu90)
            set(mpi_arguments --mpiexec "${MPIEXEC_EXECUTABLE}")
        endif()
        if(NOT variant MATCHES "_mpi$")
            target_sources(${target} PRIVATE "${app}/mpi_s.f90")
        endif()
        if(variant MATCHES "mem")
            target_sources(${target} PRIVATE "${app}/setmcp_mem.f90" "${app}/setmcp_mem_I.f90"
                "${app}/rmcdhf_mem_C.f90" "${app}/systemmem.f90" "${app}/systemmem_I.f90"
                "${app}/systemfreemem.f90" "${app}/systemfreemem_I.f90")
        endif()
        add_test(NAME descriptor_input_${variant} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_rmcdhf_descriptors.py
            --angular $<TARGET_FILE:${angular}> --loader $<TARGET_FILE:${target}> ${mpi_arguments})
        set_tests_properties(descriptor_input_${variant} PROPERTIES TIMEOUT 180)
    endif()
endforeach()

add_executable(test.descriptor_lifecycle descriptor_lifecycle.f90)
target_link_libraries_Fortran(test.descriptor_lifecycle PRIVATE 9290 mod)

if(TARGET rbiotransform AND TARGET rtransition AND TARGET rtransition_phase AND Python3_Interpreter_FOUND)
    add_test(NAME descriptor_dual_basis_pipeline COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/test_descriptor_pipeline.py --bin $<TARGET_FILE_DIR:rbiotransform>)
    set_tests_properties(descriptor_dual_basis_pipeline PROPERTIES TIMEOUT 180)
endif()
if(TARGET mpiu90)
    add_executable(test.descriptor_collectives descriptor_collectives.f90)
    target_link_libraries_Fortran(test.descriptor_collectives PRIVATE mpiu90 9290 mod)
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
        set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
        set(target test.descriptor_merge_${variant})
        if(variant MATCHES "_mpi$")
            set(suffix mpi)
        else()
            set(suffix "")
        endif()
        add_executable(${target} descriptor_merge.F90
            "${app}/mrgcsl${suffix}.f90" "${app}/mrgcsl${suffix}_I.f90"
            "${app}/ldcsl1${suffix}.f90" "${app}/ldcsl1${suffix}_I.f90"
            "${app}/ldcsl2${suffix}.f90" "${app}/ldcsl2${suffix}_I.f90"
            "${app}/merg12${suffix}.f90" "${app}/merg12${suffix}_I.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE 9290 mod)
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE DESCRIPTOR_MPI)
            target_link_libraries_Fortran(${target} PRIVATE mpiu90)
            set(mpi_arguments --mpiexec "${MPIEXEC_EXECUTABLE}")
        endif()
        add_test(NAME descriptor_merge_${variant} COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_descriptor_merge.py $<TARGET_FILE:${target}> ${mpi_arguments})
        set_tests_properties(descriptor_merge_${variant} PROPERTIES TIMEOUT 60)
    endif()
endforeach()

set(RCSFS_TEST_CLI "$ENV{RCSFS_TEST_CLI}" CACHE FILEPATH "Native rcsfs CLI used only to prepare test fixtures")
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
