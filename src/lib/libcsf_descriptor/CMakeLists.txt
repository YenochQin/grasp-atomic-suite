# Build the suite-owned Rust reader alongside the ordinary numerical libraries.
foreach(legacy RCSFS_DESCRIPTOR_SOURCE_DIR RCSFS_DESCRIPTOR_LIBRARY)
    if(DEFINED ${legacy} AND NOT "${${legacy}}" STREQUAL "")
        message(FATAL_ERROR
            "${legacy} is no longer supported. Remove it from the CMake cache; "
            "the suite builds src/lib/libcsf_descriptor into its own lib directory.")
    endif()
endforeach()
option(GRASP_CSF_READER_OFFLINE "Build the bundled reader using cached Cargo dependencies" ON)
set(GRASP_RUST_VERSION "1.99.0" CACHE STRING "Minimum Rust release version for the bundled reader")
if(NOT GRASP_RUST_VERSION MATCHES "^[0-9]+\\.[0-9]+\\.[0-9]+$")
    message(FATAL_ERROR "GRASP_RUST_VERSION must be a numeric major.minor.patch version")
endif()
find_program(GRASP_CARGO cargo REQUIRED)
find_program(GRASP_RUSTC rustc REQUIRED)
execute_process(COMMAND "${GRASP_RUSTC}" --version OUTPUT_VARIABLE rustc_version_output
    OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
string(REGEX MATCH "^rustc ([0-9]+\\.[0-9]+\\.[0-9]+)( |$)" rustc_release "${rustc_version_output}")
if(NOT rustc_release OR CMAKE_MATCH_1 VERSION_LESS GRASP_RUST_VERSION)
    message(FATAL_ERROR "Rust release ${GRASP_RUST_VERSION} or newer is required. Found ${rustc_version_output}.")
endif()
if(NOT UNIX)
    message(FATAL_ERROR "Native reader linking is currently validated on Unix; provide a platform linker implementation")
endif()

set(reader_manifest "${CMAKE_CURRENT_SOURCE_DIR}/Cargo.toml")
set(reader_target_dir "${CMAKE_CURRENT_BINARY_DIR}/cargo")
set(reader_cargo_archive "${reader_target_dir}/release/libcsf_descriptor.a")
set(reader_archive "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}/libcsf_descriptor.a")
set(cargo_network_args)
if(GRASP_CSF_READER_OFFLINE)
    list(APPEND cargo_network_args --offline)
endif()
set(reader_sources src/lib.rs src/descriptor_transport.rs src/descriptor_reader.rs src/descriptor_ffi.rs)
add_custom_command(OUTPUT "${reader_archive}"
    COMMAND "${CMAKE_COMMAND}" -E env "RUSTC=${GRASP_RUSTC}" "${GRASP_CARGO}" build
        --manifest-path "${reader_manifest}" --target-dir "${reader_target_dir}"
        --release --locked ${cargo_network_args}
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}"
    COMMAND "${CMAKE_COMMAND}" -E copy_if_different "${reader_cargo_archive}" "${reader_archive}"
    COMMAND "${CMAKE_COMMAND}" -E touch "${reader_archive}"
    DEPENDS ${reader_sources} Cargo.toml Cargo.lock
    BYPRODUCTS "${reader_cargo_archive}"
    COMMENT "Building suite CSF descriptor library"
    VERBATIM)
add_custom_target(csf_descriptor_build DEPENDS "${reader_archive}")
add_library(csf_descriptor STATIC IMPORTED GLOBAL)
set_target_properties(csf_descriptor PROPERTIES IMPORTED_LOCATION "${reader_archive}")
add_dependencies(csf_descriptor csf_descriptor_build)
find_package(Threads REQUIRED)
set_property(TARGET csf_descriptor PROPERTY INTERFACE_LINK_LIBRARIES "Threads::Threads;${CMAKE_DL_LIBS};m")
install(FILES "${reader_archive}" DESTINATION lib)
install(FILES include/rcsfs_descriptor.h DESTINATION include)
install(FILES LICENSE DESTINATION share/licenses/grasp-csf-descriptor)

get_filename_component(rustc_bin_dir "${GRASP_RUSTC}" DIRECTORY)
find_program(GRASP_RUSTDOC rustdoc HINTS "${rustc_bin_dir}" NO_DEFAULT_PATH NO_CACHE)
set(reader_test_environment "RUSTC=${GRASP_RUSTC}")
if(GRASP_RUSTDOC)
    list(APPEND reader_test_environment "RUSTDOC=${GRASP_RUSTDOC}")
endif()
add_test(NAME csf_descriptor_rust
    COMMAND "${CMAKE_COMMAND}" -E env ${reader_test_environment} "${GRASP_CARGO}" test
        --manifest-path "${reader_manifest}" --target-dir "${reader_target_dir}"
        --release --locked ${cargo_network_args} --jobs 4)
set_tests_properties(csf_descriptor_rust PROPERTIES TIMEOUT 300 PROCESSORS 4)
