! Order three active subshells and apply the tensor permutation phase.
module permuted_three_subshell_recoupling
    implicit none (type, external)

    ! Module interface
    private
    public :: rec3
contains
    subroutine rec3(first_shell_position,second_shell_position,third_shell_position,two_rank_1,two_rank_2, &
        two_rank_3,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
        ! Module dependencies
        use debug_output_settings, only: angular_debug_unit, ensure_debug_output
        use, intrinsic :: iso_fortran_env, only: real64
        use three_subshell_recoupling, only: reco3
        implicit none (type, external)

        ! Arguments
        integer, intent(in)       :: first_shell_position,second_shell_position,third_shell_position,two_rank_1, &
            two_rank_2,two_rank_3,evaluate_coefficient
        integer, intent(out)      :: coupling_allowed
        real(real64), intent(out) :: recoupling_coefficient

        ! Local declarations
        integer :: phase_sum

        if (first_shell_position == second_shell_position) then
            call ensure_debug_output()
            write(angular_debug_unit,'(5x,"ERRO IN REC")')
            stop
        end if

        ! Use the ordered three-shell routine and compensate for the tensor order.
        if((third_shell_position > first_shell_position).and.(third_shell_position > second_shell_position)) then
            if(first_shell_position-second_shell_position < 0) then
                call reco3(first_shell_position,second_shell_position,third_shell_position,two_rank_1,two_rank_2, &
                    two_rank_3,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
            else
                call reco3(second_shell_position,first_shell_position,third_shell_position,two_rank_2,two_rank_1, &
                    two_rank_3,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
                phase_sum=two_rank_1+two_rank_2-two_rank_3
                if((phase_sum/4)*4 /= phase_sum)recoupling_coefficient=-recoupling_coefficient
            end if
        else if((third_shell_position < first_shell_position).and.(third_shell_position < second_shell_position)) then
            if(first_shell_position-second_shell_position < 0) then
                call reco3(third_shell_position,first_shell_position,second_shell_position,two_rank_3,two_rank_1, &
                    two_rank_2,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
                if((two_rank_3/2)*2 /= two_rank_3)recoupling_coefficient=-recoupling_coefficient
            else
                call reco3(third_shell_position,second_shell_position,first_shell_position,two_rank_3,two_rank_2, &
                    two_rank_1,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
                phase_sum=two_rank_1+two_rank_2+two_rank_3
                if((phase_sum/4)*4 /= phase_sum)recoupling_coefficient=-recoupling_coefficient
            end if
        else
            if(first_shell_position-second_shell_position < 0)then
                call reco3(first_shell_position,third_shell_position,second_shell_position,two_rank_1,two_rank_3, &
                    two_rank_2,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
                phase_sum=two_rank_1-two_rank_2-two_rank_3
                if((phase_sum/4)*4 /= phase_sum)recoupling_coefficient=-recoupling_coefficient
            else
                call reco3(second_shell_position,third_shell_position,first_shell_position,two_rank_2,two_rank_3, &
                    two_rank_1,evaluate_coefficient,coupling_allowed,recoupling_coefficient)
                if((two_rank_1/2)*2 /= two_rank_1)recoupling_coefficient=-recoupling_coefficient
            end if
        end if
        return
    end subroutine rec3
end module permuted_three_subshell_recoupling
