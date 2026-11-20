module angular_coupling_workspace
    implicit none (type, external)

    ! Module constants
    integer, parameter :: mangm = 60
    integer, parameter :: mtriad = 12

    ! Module state
    integer :: j1(mangm)
    integer :: j2(mtriad, 3), j3(mtriad, 3)
    integer :: mja, nja
    logical :: free(mangm)
end module angular_coupling_workspace
