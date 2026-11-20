! Stores indices and values for the two angular matrix-element work lists.
module angular_matrix_element_workspace
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer                       :: inumb_int1, inumb_int2
    integer      :: int1(2000), int2(2000)
    real(real64) :: aint1(2000), aint2(2000)
end module angular_matrix_element_workspace
