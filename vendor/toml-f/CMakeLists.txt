# Suite integration for the unmodified toml-f 0.5.2 sources (MIT license).
# Keep the parser private to rwfnestimate, with no configure-time downloads.
set(srcs)
add_subdirectory(src)
add_library(suite_tomlf STATIC ${srcs})
set_target_properties(suite_tomlf PROPERTIES
    Fortran_MODULE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/suite_tomlf")
file(MAKE_DIRECTORY "${Fortran_MODULE_DIRECTORY_root}/suite_tomlf")
install(FILES LICENSE-MIT LICENSE-Apache UPSTREAM.json
    DESTINATION share/licenses/toml-f)
