add_executable(test.csfg_capacity capacity.f90)
target_link_libraries_Fortran(test.csfg_capacity PRIVATE mod_csfg mod)
add_test(NAME csfg_capacity COMMAND test.csfg_capacity)
if(Python3_Interpreter_FOUND)
    add_test(NAME csfg_file_formats
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/check.py
        $<TARGET_FILE_DIR:rcsfggenerate_csfg>)
    set_tests_properties(csfg_file_formats PROPERTIES TIMEOUT 90)
endif()
if(Python3_Interpreter_FOUND AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_test(NAME csfg_alternate_capacity
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/alternate_capacity.py
        ${CMAKE_Fortran_COMPILER})
    set_tests_properties(csfg_alternate_capacity PROPERTIES TIMEOUT 120)
endif()
