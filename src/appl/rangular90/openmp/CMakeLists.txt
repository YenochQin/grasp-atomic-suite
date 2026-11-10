# GNU -fopenmp already implies automatic/reentrant local storage. No legacy
# static-storage override or redundant -frecursive/-fautomatic is needed.
include(concurrent_sources.cmake)

# Direct objects take precedence over the legacy archives. Exactly one instance
# of each module is linked into rangular_omp; its mutable state uses OpenMP TLS.
add_library(rangular_thread_state OBJECT
    ../../../lib/libmod/m_C.f90
    ../../../lib/libmod/dumx_C.f90
    ../../../lib/libmod/trk_C.f90
    ../../../lib/libmod/buffer_C.f90
    ../../../lib/libmod/memory_man.f90)
setup_fortran_modules(rangular_thread_state)
target_link_libraries_Fortran(rangular_thread_state PRIVATE mod)
target_link_libraries(rangular_thread_state PRIVATE OpenMP::OpenMP_Fortran)

add_library(rangular_thread_kernels OBJECT ${rangular_concurrent_sources})
setup_fortran_modules(rangular_thread_kernels)
target_link_libraries_Fortran(rangular_thread_kernels PRIVATE rangular_thread_state rang90 mcp90 9290 mod)
target_include_directories(rangular_thread_kernels BEFORE PRIVATE
    "$<TARGET_PROPERTY:rangular_thread_state,Fortran_MODULE_DIRECTORY>")
target_link_libraries(rangular_thread_kernels PRIVATE OpenMP::OpenMP_Fortran)

set(rangular_omp_application_sources)
foreach(source IN LISTS rangular_common_sources)
    list(APPEND rangular_omp_application_sources "../${source}")
endforeach()
add_executable(rangular_omp rangular_parallel.f90 ${rangular_omp_application_sources})
setup_fortran_modules(rangular_omp)
target_link_libraries_Fortran(rangular_omp PRIVATE rangular_thread_kernels rangular_thread_state rang90 mcp90 9290 mod)
target_include_directories(rangular_omp BEFORE PRIVATE
    "$<TARGET_PROPERTY:rangular_thread_state,Fortran_MODULE_DIRECTORY>")
target_link_libraries(rangular_omp PRIVATE OpenMP::OpenMP_Fortran)
install(TARGETS rangular_omp DESTINATION bin/)
