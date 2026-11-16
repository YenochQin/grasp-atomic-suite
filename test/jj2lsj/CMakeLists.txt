string(REPLACE "-fno-automatic" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
add_executable(test.jj2lsj_coefficients coefficients.f90)
target_link_libraries_Fortran(test.jj2lsj_coefficients PRIVATE jj2lsj_conversion)
if(Python3_Interpreter_FOUND)
    add_test(NAME jj2lsj_parallel_regression
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check.py
        $<TARGET_FILE:jj2lsj_omp> $<TARGET_FILE:test.jj2lsj_coefficients>
        --openmp=$<BOOL:$<TARGET_PROPERTY:jj2lsj_conversion,JJ2LSJ_OPENMP>>)
    set_tests_properties(jj2lsj_parallel_regression PROPERTIES TIMEOUT 180 PROCESSORS 4)
endif()

add_executable(test.jj2lsj_ranking ranking.f90)
target_link_libraries_Fortran(test.jj2lsj_ranking PRIVATE jj2lsj_conversion)
add_test(NAME jj2lsj_stable_ranking COMMAND test.jj2lsj_ranking)
add_executable(test.jj2lsj_pointer_defaults pointer_defaults.f90)
target_link_libraries_Fortran(test.jj2lsj_pointer_defaults PRIVATE jj2lsj_conversion)
add_test(NAME jj2lsj_pointer_defaults COMMAND test.jj2lsj_pointer_defaults)
add_test(NAME jj2lsj_nonfinite_ranking COMMAND test.jj2lsj_ranking nonfinite)
set_tests_properties(jj2lsj_nonfinite_ranking PROPERTIES WILL_FAIL TRUE)

if(Python3_Interpreter_FOUND)
    add_test(NAME jj2lsj_unique_io_build COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_unique_io.py --compile
        ${CMAKE_Fortran_COMPILER} ${CMAKE_CURRENT_BINARY_DIR}/unique-io)
    set_tests_properties(jj2lsj_unique_io_build PROPERTIES FIXTURES_SETUP unique_io TIMEOUT 60)
    add_test(NAME jj2lsj_unique_io COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/check_unique_io.py ${CMAKE_CURRENT_BINARY_DIR}/unique-io/unique)
    set_tests_properties(jj2lsj_unique_io PROPERTIES FIXTURES_REQUIRED unique_io TIMEOUT 120)
endif()
