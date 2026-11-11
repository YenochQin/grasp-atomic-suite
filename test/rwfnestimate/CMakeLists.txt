if(TARGET rwfnestimate AND Python3_Interpreter_FOUND)
    foreach(element Ni Zn)
        add_test(NAME rwfnestimate_${element}_high_l
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/high_l.py
                --bindir $<TARGET_FILE_DIR:rwfnestimate> --element ${element})
        set_tests_properties(rwfnestimate_${element}_high_l PROPERTIES TIMEOUT 90)
    endforeach()
    add_test(NAME rwfnestimate_header_input
        COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/header_input.py
            --binary $<TARGET_FILE:rwfnestimate>)
    set_tests_properties(rwfnestimate_header_input PROPERTIES TIMEOUT 120 LABELS f2018_io)
endif()
