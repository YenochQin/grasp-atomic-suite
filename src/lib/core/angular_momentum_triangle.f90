! Test the triangle bounds for three angular momenta supplied as 2J+1.
module angular_momentum_triangle
    implicit none (type, external)

    ! Module interface
    private
    public :: satisfies_angular_triangle
contains
    integer function satisfies_angular_triangle (first_two_j_plus_one, second_two_j_plus_one, third_two_j_plus_one)
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: first_two_j_plus_one
        integer, intent(in) :: second_two_j_plus_one
        integer, intent(in) :: third_two_j_plus_one

        ! Local declarations
        integer :: second_minus_third
        second_minus_third = second_two_j_plus_one - third_two_j_plus_one
        if (first_two_j_plus_one>=abs(second_minus_third) + 1 .and. &
            & first_two_j_plus_one<=second_two_j_plus_one+third_two_j_plus_one-1) then
            satisfies_angular_triangle = 1
        else
            satisfies_angular_triangle = 0
        end if
        return
    end function satisfies_angular_triangle
end module angular_momentum_triangle
