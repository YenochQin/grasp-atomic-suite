module read_eol_level_weights_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_eol_level_weights
contains

    ! Interactively determines the weights for EOL calculation.
    ! It's modified to always ask the question for the weight
    !
    ! Written by Xinghong He. Last revision: 19 Mar 1999
    subroutine read_eol_level_weights(nondefault_options_flag, selected_level_count, level_weights)
        ! Module dependencies
        use calculation_file_units, only: calculation_log_unit
        use, intrinsic :: iso_fortran_env, only: real64
        use standard_io_units, only: stderr_unit, stdin_unit
        implicit none (type, external)

        ! Arguments
        integer  :: nondefault_options_flag
        integer , intent(in) :: selected_level_count
        real(real64), intent(inout) :: level_weights(selected_level_count)

        ! Local declarations
        integer :: i, weighting_mode
        real(real64) :: weight_sum

        ! Standard weights: ncmin=1, OL calculation

        if (selected_level_count == 1) then
            level_weights(1) = -1.0_real64
            return
        end if

        ! Select a method to assign level weights for ncmin > 1 case

        write (stderr_unit, *) "level weights (1 equal;  5 standard;  9 user)"

        ! Let user try 10 times to get the correct input. 10 is BIG
        ! enough since the idea here is to allow user mistakes and
        ! at the same time to avoid an infinity loop.

        do i = 1, 10
            read (stdin_unit, *) weighting_mode
            if (weighting_mode==1 .or. weighting_mode==5 .or. weighting_mode==9) exit
            write (stderr_unit, *) "Input not correct, do it again. tried=", i
        end do

        if (nondefault_options_flag==0) then
            write(calculation_log_unit,*) weighting_mode,"! level weights"
        end if

        if (i > 10) stop

        select case (weighting_mode)                         ! Equal weight
        case (1)
            level_weights = -2.0_real64
        case (5)                                   ! Standard weight
            level_weights = -1.0_real64
        case (9)                                   ! User-input weight
            do
                write (stderr_unit, *) "Enter the (relative) weights of the", selected_level_count, &
                    " levels :"
                read (stdin_unit, *) (level_weights(i),i=1,selected_level_count)

                weight_sum = 0.0_real64
                do i = 1, selected_level_count
                    if (level_weights(i) <= 0.0_real64) then
                        write (stderr_unit, *) "Weights must exceed 0;"
                        exit
                    else
                        weight_sum = weight_sum + level_weights(i)
                    end if
                end do
                if (i > selected_level_count) exit
            end do

            if (nondefault_options_flag == 0) then
                write(calculation_log_unit,*) (level_weights(i), i = 1, selected_level_count)
            end if

            weight_sum = 1.0_real64/weight_sum
            level_weights = weight_sum*level_weights

        case default
            write (stderr_unit, *) "Impossible ! Because it was guarded"
            stop
        end select

        return
    end subroutine read_eol_level_weights
end module read_eol_level_weights_interface
