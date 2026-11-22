! Parse the requested electric and magnetic multipoles and the level selection restrictions.
module read_multipole_options_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: read_multipole_options
contains

    ! Interactively  determines the  list of radiation multipolarities
    ! and parities. This is loadad into COMMON/OSC6/.
    !
    ! Written by Farid A. Parpia. Last revision: 28 Dec 1992
    subroutine read_multipole_options
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: output_unit
        use array_memory_management, only: alloc, ralloc
        use transition_off_diagonal_flags, only: m1_different_j_only_flag, e2_different_j_only_flag
        use transition_matrix_workspace, only: requested_multipole_count, requested_multipole_codes
        use yes_no_input, only: read_yes_no
        use, intrinsic :: iso_fortran_env, only: real32
        implicit none (type, external)

        ! Local declarations
        integer :: ndkp, istart, i, iend, lenth, io_status, mult
        logical :: lelec, lmagn, yes
        character :: record*256, reci
        ! Reserve space for the first requested multipole code.
        ndkp = 1
        call alloc (requested_multipole_codes, ndkp, "requested_multipole_codes", "read_multipole_options" )

        do
            write (output_unit, *) "Enter the list of transition specifications"
            write (output_unit, *) " e.g.,  E1,M2  or  E1 M2  or  E1;M2 :"
            do
                read (*, "(A)") record
                requested_multipole_count = 0
                ! Scan multipole tokens from left to right.
                istart = 0
                do i = 1, len(record)
                    reci = record(i:i)
                    if (reci/=" " .and. reci/="," .and. reci/=";") then
                        if (istart == 0) istart = i
                    else
                        if (istart /= 0) then
                            iend = i - 1
                            reci = record(istart:istart)
                            if (reci == "E") then
                                lelec = .true.
                                lmagn = .false.
                            else if (reci == "M") then
                                lelec = .false.
                                lmagn = .true.
                            else
                                write (output_unit, *) "read_multipole_options: Transitions must be of type"
                                write (output_unit, *) " E or type M; reenter ..."
                                exit
                            end if
                            lenth = iend - istart
                            if (lenth /= 1) then
                                write (output_unit, *) "read_multipole_options: Transition multipolarities"
                                write (output_unit, *) " must be integers between 1 and 9;"
                                write (output_unit, *) " reenter ..."
                                exit
                            end if
                            reci = record(iend:iend)
                            read (reci, "(1I1)", iostat=io_status) mult
                            if (io_status /= 0) then
                                write (output_unit, *) "read_multipole_options: Unable to decode multipolarity"
                                write (output_unit, *) " "//reci//"; reenter ..."
                                exit
                            end if
                            requested_multipole_count = requested_multipole_count + 1
                            if (requested_multipole_count > ndkp) then
                                call ralloc (requested_multipole_codes, requested_multipole_count, "requested_multipole_codes", &
                                    "read_multipole_options" )
                                ndkp = requested_multipole_count
                            end if
                            if (lelec) then
                                requested_multipole_codes(requested_multipole_count) = mult*(-1)**mult
                            else if (lmagn) then
                                requested_multipole_codes(requested_multipole_count) = mult*(-1)**(mult + 1)
                            end if
                            istart = 0
                        end if
                    end if
                end do
                ! An invalid token stops the scan and retries the input.
                if (i > len(record)) exit
            end do
            if (requested_multipole_count /= 0) exit
        end do
        ! Trim the multipole array to the number of accepted requests.
        if (ndkp /= requested_multipole_count) call ralloc (requested_multipole_codes, requested_multipole_count, &
            "requested_multipole_codes", "read_multipole_options")
        ! For M1 and E2, ask whether the two levels must have different J.
        do i = 1, requested_multipole_count
            if (requested_multipole_codes(i) == 1) then
                write (*, *) "M1 transitions only between levels with different J?"
                yes = read_yes_no()
                if (yes) then
                    m1_different_j_only_flag = 1
                else
                    m1_different_j_only_flag = 0
                end if
            end if
            if (requested_multipole_codes(i) /= 2) cycle
            write (*, *) "E2 transitions only between levels with different J?"
            yes = read_yes_no()
            if (yes) then
                e2_different_j_only_flag = 1
            else
                e2_different_j_only_flag = 0
            end if
        end do
        return
    end subroutine read_multipole_options
end module read_multipole_options_interface
