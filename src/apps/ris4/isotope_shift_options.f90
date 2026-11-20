module isotope_shift_options
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    !...Translated by Gediminas Gaigalas 11/18/19
    real(real64), pointer :: dens1(:), dens2(:), dens3(:),  &
        dens4(:), dens5(:), dens6(:),  &
        dens7(:), smsc1(:), smsc2(:)
end module isotope_shift_options
