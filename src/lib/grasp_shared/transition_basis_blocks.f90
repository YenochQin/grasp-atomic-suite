! Stores CSF and eigenvector counts, symmetries, and index arrays for transition bases.
module transition_basis_blocks
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ncfblk(20), nevinblk(20), ncfinblk(20), two_j(20)
    integer :: nblock
    real(real64), pointer :: idxblk(:)
    integer :: nelectot, ncftot, nwtot, nvectot, nvecsiztot, nblock1
    integer :: nblocki
    integer :: nblockf
    integer :: ncfi(10)
    integer :: ncff(10)
end module transition_basis_blocks
