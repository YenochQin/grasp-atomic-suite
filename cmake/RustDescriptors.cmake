# Explicit source or prebuilt-library dependency. Never discover a sibling checkout.
set(RCSFS_DESCRIPTOR_LIBRARY "" CACHE FILEPATH "Prebuilt native rCSFs descriptor static library")
set(RCSFS_DESCRIPTOR_SOURCE_DIR "" CACHE PATH "rCSFs source root containing native-reader/Cargo.toml")
option(RCSFS_DESCRIPTOR_OFFLINE "Use only cached/vendored Cargo dependencies" ON)
set(RCSFS_RUST_VERSION "1.99.0" CACHE STRING "Reviewed Rust compiler version for source builds")
if(RCSFS_DESCRIPTOR_SOURCE_DIR)
    find_program(RCSFS_CARGO cargo REQUIRED)
    find_program(RCSFS_RUSTC rustc REQUIRED)
    execute_process(COMMAND "${RCSFS_RUSTC}" --version OUTPUT_VARIABLE rust_version
        OUTPUT_STRIP_TRAILING_WHITESPACE COMMAND_ERROR_IS_FATAL ANY)
    string(FIND "${rust_version}" "rustc ${RCSFS_RUST_VERSION} " compiler_match)
    if(NOT compiler_match EQUAL 0)
        message(FATAL_ERROR "Expected Rust ${RCSFS_RUST_VERSION}, found ${rust_version}; select the reviewed toolchain explicitly")
    endif()
    set(reader_manifest "${RCSFS_DESCRIPTOR_SOURCE_DIR}/native-reader/Cargo.toml")
    if(NOT EXISTS "${reader_manifest}")
        message(FATAL_ERROR "RCSFS_DESCRIPTOR_SOURCE_DIR must contain native-reader/Cargo.toml")
    endif()
    if(NOT UNIX)
        message(FATAL_ERROR "Native reader linking is currently validated on Unix; provide a platform linker implementation")
    endif()
    set(reader_archive "${CMAKE_BINARY_DIR}/rust-descriptor/release/librcsfs_descriptor.a")
    set(cargo_network_args)
    if(RCSFS_DESCRIPTOR_OFFLINE)
        list(APPEND cargo_network_args --offline)
    endif()
    set(reader_sources "${RCSFS_DESCRIPTOR_SOURCE_DIR}/src/descriptor_transport.rs"
        "${RCSFS_DESCRIPTOR_SOURCE_DIR}/src/descriptor_reader.rs"
        "${RCSFS_DESCRIPTOR_SOURCE_DIR}/src/descriptor_ffi.rs")
    add_custom_command(OUTPUT "${reader_archive}"
        COMMAND "${CMAKE_COMMAND}" -E env "RUSTC=${RCSFS_RUSTC}" "${RCSFS_CARGO}" build --manifest-path "${reader_manifest}"
            --target-dir "${CMAKE_BINARY_DIR}/rust-descriptor" --release --locked ${cargo_network_args}
        DEPENDS ${reader_sources} "${reader_manifest}"
            "${RCSFS_DESCRIPTOR_SOURCE_DIR}/native-reader/lib.rs"
            "${RCSFS_DESCRIPTOR_SOURCE_DIR}/native-reader/Cargo.lock"
        COMMENT "Building native descriptor reader (no Python, CLI or CSF generator)"
        VERBATIM)
    add_custom_target(rcsfs_descriptor_build DEPENDS "${reader_archive}")
elseif(EXISTS "${RCSFS_DESCRIPTOR_LIBRARY}")
    set(reader_archive "${RCSFS_DESCRIPTOR_LIBRARY}")
else()
    message(FATAL_ERROR "Set RCSFS_DESCRIPTOR_SOURCE_DIR or RCSFS_DESCRIPTOR_LIBRARY; native descriptors are required")
endif()
add_library(rcsfs_descriptor_native STATIC IMPORTED GLOBAL)
set_target_properties(rcsfs_descriptor_native PROPERTIES IMPORTED_LOCATION "${reader_archive}")
if(TARGET rcsfs_descriptor_build)
    add_dependencies(rcsfs_descriptor_native rcsfs_descriptor_build)
endif()
find_package(Threads REQUIRED)
set_property(TARGET rcsfs_descriptor_native PROPERTY INTERFACE_LINK_LIBRARIES "Threads::Threads;${CMAKE_DL_LIBS};m")
install(FILES "${reader_archive}" DESTINATION lib RENAME librcsfs_descriptor.a)
