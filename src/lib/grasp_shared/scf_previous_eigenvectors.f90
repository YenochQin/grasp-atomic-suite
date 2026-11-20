module scf_previous_eigenvectors
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    real(real64), pointer :: previous_csf_mixing_coefficients(:)
end module scf_previous_eigenvectors
