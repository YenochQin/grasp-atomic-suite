module rk_fixture
    implicit none (type, external)
    integer :: calls=0
end module rk_fixture

module slater_radial_integral
    implicit none (type, external)
    private
    public :: calculate_slater_radial_integral
contains

function calculate_slater_radial_integral(ia,ib,ic,id,k) result(value)
    use, intrinsic :: iso_fortran_env, only: real64
    use rk_fixture, only: calls
    implicit none (type, external)
    integer,intent(in) :: ia,ib,ic,id,k
    real(real64) :: value
    if(min(ia,ib,ic,id)<1.or.k<0) error stop "Invalid Slater fixture indices"
    calls=calls+1
    value=1.0_real64
end function calculate_slater_radial_integral
end module slater_radial_integral
