! Evaluate a Racah coefficient by a descending finite-sum recurrence.
module racah_coefficients
    implicit none (type, external)

    ! Module interface
    private
    public :: dracah
contains
    subroutine dracah(two_j_1, two_j_2, two_j_3, two_j_4, two_j_5, two_j_6, racah_coefficient)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use log_factorial_cache, only: log_factorials => gam, factorial_table_capacity => mfact
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: two_j_1
        integer , intent(in) :: two_j_2
        integer , intent(in) :: two_j_3
        integer , intent(in) :: two_j_4
        integer , intent(in) :: two_j_5
        integer , intent(in) :: two_j_6
        real(real64) , intent(out) :: racah_coefficient

        ! Local declarations
        integer :: triangle_sum_1, triangle_sum_2, triangle_sum_3, triangle_sum_4, quadrangle_sum_1, &
            quadrangle_sum_2, quadrangle_sum_3, sum_start, sum_end, terms_evaluated, term_iteration, &
            descending_term_index

        ! Form doubled triangle sums for selection, then use their halves in the Racah sum.
        triangle_sum_1 = two_j_1 + two_j_2 + two_j_5
        triangle_sum_2 = two_j_3 + two_j_4 + two_j_5
        triangle_sum_3 = two_j_1 + two_j_3 + two_j_6
        triangle_sum_4 = two_j_2 + two_j_4 + two_j_6
        racah_coefficient = 0.0_real64
        if (2*max(max(two_j_1,two_j_2),two_j_5) - triangle_sum_1>0 .or. mod(triangle_sum_1,2)/=0) return
        if (2*max(max(two_j_3,two_j_4),two_j_5) - triangle_sum_2>0 .or. mod(triangle_sum_2,2)/=0) return
        if (2*max(max(two_j_1,two_j_3),two_j_6) - triangle_sum_3>0 .or. mod(triangle_sum_3,2)/=0) return
        if (2*max(max(two_j_2,two_j_4),two_j_6) - triangle_sum_4>0 .or. mod(triangle_sum_4,2)/=0) return
        triangle_sum_1 = triangle_sum_1/2
        triangle_sum_2 = triangle_sum_2/2
        triangle_sum_3 = triangle_sum_3/2
        triangle_sum_4 = triangle_sum_4/2
        quadrangle_sum_1 = (two_j_1 + two_j_2 + two_j_3 + two_j_4)/2
        quadrangle_sum_2 = (two_j_1 + two_j_4 + two_j_5 + two_j_6)/2
        quadrangle_sum_3 = (two_j_2 + two_j_3 + two_j_5 + two_j_6)/2

        ! The largest triangle sum and smallest quadrangle sum bound the summation.
        sum_start = max(max(max(triangle_sum_1,triangle_sum_2),triangle_sum_3),triangle_sum_4) + 1
        sum_end = min(min(quadrangle_sum_1,quadrangle_sum_2),quadrangle_sum_3) + 1
        racah_coefficient = 1.0_real64
        terms_evaluated = 0
        if (sum_start /= sum_end) then
            sum_start = sum_start + 1

            ! Evaluate the alternating sum by a descending recurrence in adjacent-term ratios.
            do term_iteration = sum_start, sum_end
                descending_term_index = sum_end - terms_evaluated
                racah_coefficient = 1.0_real64 - racah_coefficient*dble(descending_term_index*(quadrangle_sum_1 - &
                    descending_term_index + 2)*(quadrangle_sum_2 - descending_term_index + 2)*(quadrangle_sum_3 - &
                    descending_term_index + 2&
                    ))/dble((descending_term_index - 1 - triangle_sum_1)*(descending_term_index - 1 - triangle_sum_2)* &
                    (descending_term_index - 1 - triangle_sum_3)*(descending_term_index - 1 - triangle_sum_4))
                terms_evaluated = terms_evaluated + 1
            end do
            sum_start = sum_start - 1
        end if
        racah_coefficient = racah_coefficient*(-1.0_real64)**(quadrangle_sum_1 + sum_start + 1)* &
            exp((log_factorials(sum_start+1)-log_factorials(sum_start-triangle_sum_1)-log_factorials&
            (sum_start-triangle_sum_2)-log_factorials(sum_start-triangle_sum_3)-log_factorials(sum_start- &
            triangle_sum_4)-log_factorials(quadrangle_sum_1+2-sum_start)-log_factorials(quadrangle_sum_2+2- &
            sum_start)&
            -log_factorials(quadrangle_sum_3+2-sum_start))+(log_factorials(triangle_sum_1+1-two_j_1)+ &
            log_factorials(triangle_sum_1+1-two_j_2)+log_factorials(triangle_sum_1+1-two_j_5)- &
            log_factorials(triangle_sum_1+2)+log_factorials(&
            triangle_sum_2+1-two_j_3)+log_factorials(triangle_sum_2+1-two_j_4)+log_factorials(triangle_sum_2+1-two_j_5)- &
            log_factorials(triangle_sum_2+2)+log_factorials(triangle_sum_3+1-two_j_1)+ &
            log_factorials(triangle_sum_3+1-two_j_3)+log_factorials(&
            triangle_sum_3+1-two_j_6)-log_factorials(triangle_sum_3+2)+log_factorials(triangle_sum_4+1-two_j_2)+ &
            log_factorials(triangle_sum_4+1-two_j_4)+log_factorials(triangle_sum_4+1-two_j_6)- &
            log_factorials(triangle_sum_4+2))*&
            0.5_real64)
        return
    end subroutine dracah
end module racah_coefficients
