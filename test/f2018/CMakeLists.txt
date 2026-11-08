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
