add_executable(test.f2018_davidson davidson_contracts.F90)
target_link_libraries_Fortran(test.f2018_davidson PRIVATE dvd90 9290 mod)
target_link_libraries(test.f2018_davidson PRIVATE ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
add_test(NAME f2018_davidson_contracts COMMAND test.f2018_davidson)
add_test(NAME f2018_davidson_early COMMAND test.f2018_davidson early)
set_tests_properties(f2018_davidson_contracts f2018_davidson_early PROPERTIES TIMEOUT 30)
foreach(variant rci rci_mpi)
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
    foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
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
foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
    if(variant MATCHES "_mpi$" AND NOT TARGET mpiu90)
        continue()
    endif()
    set(directory ${variant})
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
    set(angular_app "${PROJECT_SOURCE_DIR}/src/appl/rangular_mpi")
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

# Input failure tests call the production routines with valid, empty and
# truncated records; a single damaged MPI rank must abort without a barrier.
if(Python3_Interpreter_FOUND)
    add_test(NAME f2018_input_precision
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_input_precision.py
        ${CMAKE_Fortran_COMPILER} "${Fortran_MODULE_DIRECTORY_root}/mod")
    set_tests_properties(f2018_input_precision PROPERTIES TIMEOUT 60)
    foreach(variant rmcdhf rmcdhf_mem rmcdhf_mpi rmcdhf_mem_mpi)
        if(NOT TARGET ${variant})
            continue()
        endif()
        set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
        set(target test.mcp_io_${variant})
        add_executable(${target} mcp_reader_io.F90 "${app}/setmcp.f90" "${app}/setmcp_I.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE mod 9290)
        target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace -finit-integer=0)
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE MCP_MPI)
            target_link_libraries_Fortran(${target} PRIVATE mpiu90)
            set(mpi_arguments --mpiexec "${MPIEXEC_EXECUTABLE}")
        else()
            target_sources(${target} PRIVATE "${app}/mpi_s.f90")
        endif()
        if(variant MATCHES "mem")
            target_compile_definitions(${target} PRIVATE MCP_MEM)
            target_sources(${target} PRIVATE "${app}/setmcp_mem.f90" "${app}/setmcp_mem_I.f90"
                "${app}/rmcdhf_mem_C.f90" "${app}/systemmem_I.f90" "${app}/systemfreemem_I.f90")
        else()
            target_sources(${target} PRIVATE "${app}/setcof.f90" "${app}/setcof_I.f90"
                "${app}/dsubrs_I.f90" "${app}/fco_I.f90" "${app}/gco_I.f90")
        endif()
        add_test(NAME f2018_mcp_io_${variant}
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mcp_io.py
            $<TARGET_FILE:${target}> --reader ${mpi_arguments})
        set_tests_properties(f2018_mcp_io_${variant} PROPERTIES TIMEOUT 240)
        if(variant MATCHES "_mpi$")
            set_tests_properties(f2018_mcp_io_${variant} PROPERTIES PROCESSORS 2)
        endif()
        if(variant MATCHES "mem")
            set(target test.mcp_memory_${variant})
            add_executable(${target} mcp_memory_io.F90 "${app}/setmcp_mem.f90" "${app}/setmcp_mem_I.f90"
                "${app}/rmcdhf_mem_C.f90" "${app}/systemmem_I.f90" "${app}/systemfreemem_I.f90")
            setup_fortran_modules(${target})
            target_link_libraries_Fortran(${target} PRIVATE mod)
            target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace -finit-integer=0)
            if(variant MATCHES "_mpi$")
                target_compile_definitions(${target} PRIVATE MCP_MPI)
                target_link_libraries_Fortran(${target} PRIVATE mpiu90)
            else()
                target_sources(${target} PRIVATE "${app}/mpi_s.f90")
            endif()
            add_test(NAME f2018_mcp_memory_${variant}
                COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_mcp_io.py
                $<TARGET_FILE:${target}> ${mpi_arguments})
            set_tests_properties(f2018_mcp_memory_${variant} PROPERTIES TIMEOUT 180)
            if(variant MATCHES "_mpi$")
                set_tests_properties(f2018_mcp_memory_${variant} PROPERTIES PROCESSORS 2)
            endif()
        endif()
    endforeach()
    foreach(variant rbiotransform rbiotransform_mpi)
        if(NOT TARGET ${variant})
            continue()
        endif()
        set(app "${PROJECT_SOURCE_DIR}/src/appl/${variant}")
        set(target test.csf_record_${variant})
        add_executable(${target} csf_record_io.F90 "${app}/lodcslBio.f90" "${app}/lodcslBio_I.f90")
        setup_fortran_modules(${target})
        target_link_libraries_Fortran(${target} PRIVATE mod 9290)
        target_compile_options(${target} PRIVATE -fcheck=all -fbacktrace)
        set(mpi_arguments)
        if(variant MATCHES "_mpi$")
            target_compile_definitions(${target} PRIVATE MCP_MPI)
            target_link_libraries_Fortran(${target} PRIVATE mpiu90)
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
    if(TARGET mpiu90)
        add_executable(test.cpath_io cpath_io.f90)
        target_link_libraries_Fortran(test.cpath_io PRIVATE mpiu90 mod)
        add_test(NAME f2018_cpath_io_mpi
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check_cpath_io.py
            "${MPIEXEC_EXECUTABLE}" $<TARGET_FILE:test.cpath_io>)
        set_tests_properties(f2018_cpath_io_mpi PROPERTIES TIMEOUT 60 PROCESSORS 2)
    endif()
endif()
