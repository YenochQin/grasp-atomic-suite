module dpbdt_i
    implicit none (type, external)
    private
    public :: dpbdt
contains
subroutine DPBDT(j)
    use grid_c, only: N => radial_point_count
    use tatb_c, only: TA, TB
    implicit none (type, external)
    integer :: j
    TA(:N) = 0
    TB(:N) = 0
end subroutine
end module dpbdt_i

module quad_i
    implicit none (type, external)
    private
    public :: quad
contains

subroutine QUAD(result)
    use, intrinsic :: iso_fortran_env, only: real64
    use tatb_c, only: TA, MTP
    implicit none (type, external)
    real(real64), intent(out) :: result
    integer :: i
    result = 0
    do i = 2, MTP
        result = result+real(i,real64)*TA(i)
    end do
end subroutine
end module quad_i
