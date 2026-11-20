# Historical readers are test oracles only; never linked into production applications.
add_library(csf_text_oracle STATIC cslh.f90 lodcsh.f90 lodcsh2.f90 setcsll.f90
    prsrsl.f90 prsrcn.f90 parsjl.f90)
setup_fortran_modules(csf_text_oracle)
target_link_libraries_Fortran(csf_text_oracle PRIVATE core grasp_shared)
if(TARGET rangular)
    add_executable(csf_descriptor_verify csf_descriptor_verify.f90)
    setup_fortran_modules(csf_descriptor_verify)
    set_target_properties(csf_descriptor_verify PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}/bin")
    target_link_libraries_Fortran(csf_descriptor_verify PRIVATE csf_text_oracle core grasp_shared)
endif()
