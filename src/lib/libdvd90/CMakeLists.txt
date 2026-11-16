add_library(dvd90 STATIC
    adds.f90
    dvdrvr.f90
    initdvd.f90
    mgs_nrm.f90
    multbc.f90
    newvec.f90
    ovflow.f90
    tstsel.f90
    dvdson.f90
    gdvd.f90
    iniest.f90
)
setup_fortran_modules(dvd90)
target_link_libraries_Fortran(dvd90 PRIVATE 9290 mod)
install(TARGETS dvd90 DESTINATION lib/)
