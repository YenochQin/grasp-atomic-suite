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
    if(TARGET rangular_mpi AND MPIEXEC_EXECUTABLE)
        add_test(NAME rangular_mpi_descriptors
            COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rangular/test_mpi_descriptors.py"
            --rangular $<TARGET_FILE:rangular_mpi> --mpiexec "${MPIEXEC_EXECUTABLE}"
            --numproc-flag=${MPIEXEC_NUMPROC_FLAG})
        set_tests_properties(rangular_mpi_descriptors PROPERTIES TIMEOUT 180 PROCESSORS 4
            LABELS "f2018_io;f2018_mpi")
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
        if(TARGET rangular_omp)
            add_test(NAME rangular_openmp_regression
                COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rangular/test_openmp.py"
                --serial $<TARGET_FILE:rangular> --parallel $<TARGET_FILE:rangular_omp>)
            set_tests_properties(rangular_openmp_regression PROPERTIES TIMEOUT 180 PROCESSORS 8
                LABELS "f2018_io;f2018_serial")
            add_test(NAME rangular_openmp_closure
                COMMAND "${Python3_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/rangular/audit_thread_closure.py"
                "${PROJECT_BINARY_DIR}")
            set_tests_properties(rangular_openmp_closure PROPERTIES TIMEOUT 30)
        endif()
    endif()
endif()

# Independent synthetic NP tests; no GRASP preparation tools or atomic data.
add_subdirectory(npol)
add_subdirectory(jj2lsj)
add_subdirectory(standard_grasp)
add_subdirectory(rwfnestimate)
add_subdirectory(f2018)

if(TARGET mod_csfg)
    add_subdirectory(graspg)
endif()
