! Sum the alternating inverse-square and inverse-cube series used by the Fermi nuclear potential.
module alternating_exponential_series
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_alternating_exponential_series
contains
    subroutine evaluate_alternating_exponential_series(exponent_argument, inverse_square_sum, inverse_cube_sum)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        real(real64), intent(in) :: exponent_argument
        real(real64), intent(out) :: inverse_square_sum
        real(real64), intent(out) :: inverse_cube_sum

        ! Local declarations
        integer :: term_index
        real(real64) :: alternating_sign, term_index_real, inverse_term_index, exponential_weight, inverse_square_term, &
            & inverse_cube_term, previous_square_sum
        term_index = 0
        inverse_square_sum = 0.0_real64
        inverse_cube_sum = 0.0_real64
        alternating_sign = 1.0_real64
        term_index = term_index + 1
        term_index_real = dble(term_index)
        inverse_term_index = 1.0_real64/term_index_real
        alternating_sign = -alternating_sign
        exponential_weight = exp(term_index_real*exponent_argument)
        inverse_square_term = alternating_sign*exponential_weight*inverse_term_index*inverse_term_index
        inverse_cube_term = inverse_square_term*inverse_term_index
        previous_square_sum = inverse_square_sum
        inverse_square_sum = inverse_square_sum + inverse_square_term
        inverse_cube_sum = inverse_cube_sum + inverse_cube_term
        ! Stop when another term no longer changes the magnitude of the square sum in working precision.
        do while(abs(inverse_square_sum) /= abs(previous_square_sum))
            term_index = term_index + 1
            term_index_real = dble(term_index)
            inverse_term_index = 1.0_real64/term_index_real
            alternating_sign = -alternating_sign
            exponential_weight = exp(term_index_real*exponent_argument)
            inverse_square_term = alternating_sign*exponential_weight*inverse_term_index*inverse_term_index
            inverse_cube_term = inverse_square_term*inverse_term_index
            previous_square_sum = inverse_square_sum
            inverse_square_sum = inverse_square_sum + inverse_square_term
            inverse_cube_sum = inverse_cube_sum + inverse_cube_term
        end do
        return
    end subroutine evaluate_alternating_exponential_series
end module alternating_exponential_series
