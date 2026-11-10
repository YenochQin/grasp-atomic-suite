# Isolated, opt-in GRASPG/CSFG extension. It is outside the active ordinary-GRASP
# and rangular OpenMP work. Reuse common libraries/parameters; do not copy them.
add_subdirectory("src/lib/libmod_csfg")
add_subdirectory("src/lib/lib9290_csfg")
add_subdirectory("src/lib/librang90_csfg")
foreach(app jj2lsj90_csfg rwfnestimate90_csfg rmixaccumulate_csfg
            rcsfggenerate90_csfg rcsfginteract90_csfg rcsfgzerofirst90_csfg)
    add_subdirectory("src/appl/${app}")
endforeach()
add_subdirectory("src/tool/graspg")
if(TARGET mpiu90)
    add_subdirectory("src/lib/libdvd90mpi")
    add_subdirectory("src/lib/libme_csfg")
    foreach(app rangular90_csfg_mpi rci90_csfg_mpi
                rci90_block_csfg_mpi rmcdhf90_csfg_mpi)
        add_subdirectory("src/appl/${app}")
    endforeach()
endif()
