module read_al_level_weights_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_al_level_weights
contains

    !   Interactively determines the weights.
    !
    ! Written by Xinghong He. Last revision: 19 Mar 1999
    subroutine read_al_level_weights(csf_count, level_weights)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use standard_io_units, only: stderr_unit, stdin_unit

        ! Procedure interfaces
        use basis_csf_two_j_plus_one_interface, only: basis_csf_two_j_plus_one
        implicit none (type, external)

        ! Arguments
        integer , intent(in) :: csf_count
        real(real64), intent(inout) :: level_weights(csf_count)

        ! Local declarations
        integer :: i, weighting_mode
        real(real64) :: weight_sum, statistical_weight

        ! Select a method to assign level weights for NCF > 1 case

        write (stderr_unit, *) "level weights (1 equal;  5 standard;  9 user)"

        ! Let user try 10 times to get the correct input. 10 is BIG
        ! enough since the idea here is to allow user mistakes and
        ! at the same time to avoid an infinity loop.

        do i = 1, 10
            read (stdin_unit, *) weighting_mode
            if (weighting_mode==1 .or. weighting_mode==5 .or. weighting_mode==9) exit
            write (stderr_unit, *) "Input not correct, do it again. tried=", i
        end do

        if (i > 10) error stop "Must be running un-attended"

        select case (weighting_mode)                         ! Equal weight
        case (1)
            level_weights = 1.0_real64
            weight_sum = dble(csf_count)
        case (5)                                   ! Standard weight
            weight_sum = 0.0_real64
            do i = 1, csf_count
                statistical_weight = dble(basis_csf_two_j_plus_one(i))
                level_weights(i) = statistical_weight
                weight_sum = weight_sum + statistical_weight
            end do
        case (9)                                   ! User-input weight
            do
                write (stderr_unit, *) "Enter the (relative) weights of the", csf_count, &
                    " levels :"
                read (stdin_unit, *) (level_weights(i),i=1,csf_count)

                weight_sum = 0.0_real64
                do i = 1, csf_count
                    if (level_weights(i) <= 0.0_real64) then
                        write (stderr_unit, *) "Weights must exceed 0;"
                        exit
                    else
                        weight_sum = weight_sum + level_weights(i)
                    end if
                end do
                if (i > csf_count) exit
            end do

        case default
            write (stderr_unit, *) "Impossible ! Because it was guarded"
            error stop
        end select

        weight_sum = 1.0_real64/weight_sum
        level_weights = weight_sum*level_weights

        return
    end subroutine read_al_level_weights
end module read_al_level_weights_interface
