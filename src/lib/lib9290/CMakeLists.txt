add_library(9290 STATIC
    csf_rust_ffi.f90
    rmcdhf_descriptor_input.f90
    csf_basis_input.f90
    csf_basis_merge.f90
    csf_block_input.f90
    csf_descriptor_input.f90
    suite_system.f90
    npol_input.f90
    npol_kernel.f90
    npol_init.f90
    alcbuf.f90
    arctan.f90
    calen.f90
    cgamma.f90
    clrx.f90
    convrt.f90
    convrt2.f90
    convrt_double.f90
    speak.f90
    cord.f90
    count.f90
    cre.f90
    openfl.f90
    dcbsrw.f90
    dinit.f90
    dmerge.f90
    dpbdt.f90
    draw.f90
    es.f90
    factt.f90
    ldigit.f90
    getrsl.f90
    getyn.f90
    iq.f90
    ichkq1.f90
    ichkq2.f90
    ichop.f90
    icopy.f90
    iniest2.f90
    interp.f90
    quad.f90
    rint.f90
    intrpq.f90
    ispar.f90
    items.f90
    itjpo.f90
    itrig.f90
    jcup.f90
    jqs.f90
    pack.f90
    lodiso.f90
    orthsc.f90
    lodrwf.f90
    lodstate.f90
    ltab.f90
    nucpot.f90
    posfile.f90
    radgrd.f90
    rinti.f90
    setcon.f90
    setcsla.f90
    setiso.f90
    setj.f90
    setmc.f90
    setpot.f90
    setqic.f90
    setqna.f90
    setrwfa.f90
    skrc.f90
    yzk.f90
    slater.f90
    spicmv2.f90
    start.f90
    starttime.f90
    stoptime.f90
)
setup_fortran_modules(9290)
target_link_libraries(9290 PUBLIC csf_descriptor)
target_link_libraries_Fortran(9290 PRIVATE mod)
target_link_libraries(9290 PRIVATE ${BLAS_LIBRARIES} ${BLAS_LINKER_FLAGS})
target_link_libraries(9290 PRIVATE ${LAPACK_LIBRARIES} ${LAPACK_LINKER_FLAGS})
install(TARGETS 9290 DESTINATION lib/)
# Check NP production routines as well as their test callers in Debug builds.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set_source_files_properties(npol_input.f90 npol_kernel.f90 npol_init.f90 PROPERTIES
        COMPILE_OPTIONS "$<$<CONFIG:Debug>:-fcheck=all;-fbacktrace>")
endif()
