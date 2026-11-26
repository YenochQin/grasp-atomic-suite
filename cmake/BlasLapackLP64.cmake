# Suite callers pass default Fortran integers to the BLAS/LAPACK interfaces.
if(DEFINED BLA_SIZEOF_INTEGER AND NOT BLA_SIZEOF_INTEGER MATCHES "^(4|ANY|)$")
    message(FATAL_ERROR "Suite requires LP64 BLAS/LAPACK (BLA_SIZEOF_INTEGER=4, 32-bit integers). "
        "Requested BLA_SIZEOF_INTEGER=${BLA_SIZEOF_INTEGER}.")
endif()
set(BLA_SIZEOF_INTEGER 4)

# The CMake Find modules use implicit, argument-free calls to locate symbols.
# Relax only those searches and retain the configured flags for the ABI probes.
block(PROPAGATE BLAS_LIBRARIES BLAS_LINKER_FLAGS LAPACK_LIBRARIES LAPACK_LINKER_FLAGS)
    string(REPLACE "-Werror=implicit-interface" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    string(REPLACE "-Werror=implicit-procedure" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endblock()

if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    message(FATAL_ERROR "The Suite LP64 BLAS/LAPACK ABI check must execute target code. "
        "Set CMAKE_CROSSCOMPILING_EMULATOR for a cross build.")
endif()

function(grasp_verify_lp64 component)
    set(probe_dir "${CMAKE_BINARY_DIR}/CMakeFiles/grasp-linalg")
    set(probe "${probe_dir}/${component}${CMAKE_EXECUTABLE_SUFFIX}")
    file(MAKE_DIRECTORY "${probe_dir}")
    set(CMAKE_TRY_COMPILE_TARGET_TYPE EXECUTABLE)
    set(CMAKE_TRY_COMPILE_CONFIGURATION "${CMAKE_BUILD_TYPE}")
    try_compile(linked
        SOURCES
            "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/../src/lib/grasp_shared/blas_lapack_interfaces.f90"
            "${CMAKE_CURRENT_FUNCTION_LIST_DIR}/blas_lapack_lp64_probe.F90"
        COMPILE_DEFINITIONS "-DGRASP_${component}_PROBE"
        LINK_LIBRARIES ${ARGN}
        COPY_FILE "${probe}"
        LOG_DESCRIPTION "Link the Suite ${component} LP64 ABI probe"
        OUTPUT_VARIABLE link_output
        NO_CACHE)
    set(log "${probe_dir}/${component}.log")
    if(NOT linked)
        file(WRITE "${log}" "${link_output}")
        message(FATAL_ERROR "Suite ${component} LP64 ABI check could not link the selected libraries. See ${log}")
    endif()
    execute_process(
        COMMAND "${CMAKE_COMMAND}" -E env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
            ${CMAKE_CROSSCOMPILING_EMULATOR} "${probe}"
        RESULT_VARIABLE status OUTPUT_VARIABLE output ERROR_VARIABLE diagnostics
        TIMEOUT 15)
    file(WRITE "${log}" "${link_output}\nExit: ${status}\n${output}\n${diagnostics}")
    # Some implementations stop successfully after XERBLA. Require the probe's completion marker.
    if(NOT status STREQUAL "0" OR NOT output MATCHES "Suite LP64 ABI verified")
        string(STRIP "${output}\n${diagnostics}" report)
        message(FATAL_ERROR "Suite ${component} LP64 ABI check rejected the selected libraries "
            "(probe exit: ${status}). Select libraries with 32-bit Fortran integer arguments.\n${report}\nSee ${log}")
    endif()
    message(STATUS "Suite ${component} LP64 ABI verified")
endfunction()

grasp_verify_lp64(BLAS BLAS::BLAS)
# FindLAPACK's target already carries BLAS. Use the same ordering as grasp_shared.
grasp_verify_lp64(LAPACK LAPACK::LAPACK)
