module mcp_coefficient_storage
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    integer :: ncoeff, nintg
    ! intgrl holds packed orbital pairs. intptr gives the last coefficient for each pair.
    ! jann and jbnn identify the two CSFs for each entry of cnn.
    integer, pointer :: jann(:), jbnn(:), intgrl(:), intptr(:)
    real(real64), pointer :: cnn(:)
end module mcp_coefficient_storage
