if(Python3_Interpreter_FOUND AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_test(NAME standard_grasp_alternate_capacity
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/alternate_capacity.py
        ${CMAKE_Fortran_COMPILER})
    set_tests_properties(standard_grasp_alternate_capacity PROPERTIES TIMEOUT 120)
endif()
