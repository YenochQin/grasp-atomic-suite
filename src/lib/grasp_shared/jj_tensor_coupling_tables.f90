! Stores the bra and ket subshell quantum numbers used by the spin-angular tensor routines.
module jj_tensor_coupling_tables
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! perko1 defines the slot layouts. The d arrays hold ket data and the k arrays hold bra data.
    integer      :: id1(7),id2(7),ik1(7),ik2(7),id3(7),id4(7),ik3(7),ik4(7)
    real(real64) :: bd1(3),bd2(3),bk1(3),bk2(3),bd3(3),bd4(3),bk3(3),bk4(3)
end module jj_tensor_coupling_tables
