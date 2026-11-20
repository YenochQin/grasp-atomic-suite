! Stores the selected level counts and indices for the active, final, and initial bases.
module selected_level_indices
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: nvec, nvecmx
    integer :: nvecff, nvecmxff
    integer :: nvecii, nvecmxii
    real(real64) :: pnivecii
    integer, pointer :: ivec(:), ivecff(:), ivecii(:)
end module selected_level_indices
