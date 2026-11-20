! Evaluate the relativistic 3-j factor for two kappa values and a tensor rank using log-factorials.
module relativistic_three_j
    implicit none (type, external)

    ! Module interface
    private
    public :: relativistic_three_j_value
contains
    real(real64) function relativistic_three_j_value (kappa_a, tensor_rank, kappa_b)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use log_factorial_cache, only: log_factorials
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: kappa_a, tensor_rank, kappa_b

        ! Local declarations
        integer :: abs_kappa_a, abs_kappa_b, abs_kappa_sum, abs_kappa_sum_plus_rank, abs_kappa_a_minus_b, abs_kappa_b_minus_a
        real(real64) :: log_magnitude
        abs_kappa_a = abs(kappa_a)
        abs_kappa_b = abs(kappa_b)
        ! The triangle bounds use j=abs(kappa)-1/2.
        if (abs(abs_kappa_a - abs_kappa_b)<=tensor_rank .and. abs_kappa_a+abs_kappa_b-1>=tensor_rank) then
            log_magnitude = -log(dble(abs_kappa_a*abs_kappa_b))
            abs_kappa_sum = abs_kappa_a + abs_kappa_b
            abs_kappa_sum_plus_rank = abs_kappa_sum + tensor_rank
            abs_kappa_a_minus_b = abs_kappa_a - abs_kappa_b
            abs_kappa_b_minus_a = abs_kappa_b - abs_kappa_a
            log_magnitude = 0.5_real64*(log_magnitude + &
                & log_factorials(abs_kappa_sum-tensor_rank)+log_factorials(abs_kappa_a_minus_b+tensor_rank+1)+ &
                & log_factorials(abs_kappa_b_minus_a+tensor_rank+1)-&
                log_factorials(abs_kappa_sum_plus_rank+1))
            if (mod(abs_kappa_sum_plus_rank,2) == 0) then
                if (mod(3*abs_kappa_sum_plus_rank/2,2) == 0) then
                    relativistic_three_j_value = 1.0_real64
                else
                    relativistic_three_j_value = -1.0_real64
                end if
                log_magnitude = log_magnitude + log_factorials((abs_kappa_sum_plus_rank+2)/2) - &
                    & log_factorials((abs_kappa_sum-tensor_rank)/2) - log_factorials((abs_kappa_a_minus_b+&
                    tensor_rank+2)/2) - log_factorials((abs_kappa_b_minus_a+tensor_rank+2)/2)
            else
                if (mod((3*abs_kappa_sum_plus_rank - 1)/2,2) == 0) then
                    relativistic_three_j_value = 1.0_real64
                else
                    relativistic_three_j_value = -1.0_real64
                end if
                log_magnitude = log_magnitude + log_factorials((abs_kappa_sum_plus_rank+1)/2) - &
                    & log_factorials((abs_kappa_sum-tensor_rank+1)/2) - log_factorials((&
                    abs_kappa_a_minus_b+tensor_rank+1)/2) - log_factorials((abs_kappa_b_minus_a+tensor_rank+1)/2)
            end if
            relativistic_three_j_value = relativistic_three_j_value*exp(log_magnitude)
        else
            relativistic_three_j_value = 0.0_real64
        end if
        return
    end function relativistic_three_j_value
end module relativistic_three_j
