!***********************************************************************
!                                                                      *
subroutine read_multipole_options
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !   Interactively  determines the  list of radiation multipolarities   *
    !   and parities. This is loadad into COMMON/OSC6/.                    *
    !                                                                      *
    !   Call(s) to: [LIB92]: ALLOC, RALLOC.                                *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    use array_memory_management, only: alloc, ralloc
    use transition_off_diagonal_flags, only: noffd1, noffd2
    use transition_matrix_workspace, only: nkp, kp
    use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_integer, mpi_rank
    use yes_no_input, only: getyn
    use, intrinsic :: iso_fortran_env, only: real32
    implicit none (type, external)

    ! Local declarations
    integer :: ndkp, istart, i, iend, lenth, io_status, mult
    real(real32) :: number_text(3)
    logical :: lelec, lmagn, yes
    character :: record*256, reci
    !-----------------------------------------------
    !
    !
    !   Initial allocation for PNTRKP
    !
    ndkp = 1
    call alloc (kp, ndkp, "KP", "GETRMP" )
    if (mpi_rank == 0) then
        !
        !   Entry message
        !
        do
            write (output_unit, *) "Enter the list of transition specifications"
            write (output_unit, *) " e.g.,  E1,M2  or  E1 M2  or  E1;M2 :"
            do
                read (*, "(A)") record
                nkp = 0
                !
                !   Parse RECORD from left to right
                !
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
                                write (output_unit, *) "GETRMP: Transitions must be of type"
                                write (output_unit, *) " E or type M; reenter ..."
                                exit
                            end if
                            lenth = iend - istart
                            if (lenth /= 1) then
                                write (output_unit, *) "GETRMP: Transition multipolarities"
                                write (output_unit, *) " must be integers between 1 and 9;"
                                write (output_unit, *) " reenter ..."
                                exit
                            end if
                            reci = record(iend:iend)
                            read (reci, "(1I1)", iostat=io_status) mult
                            if (io_status /= 0) then
                                write (output_unit, *) "GETRMP: Unable to decode multipolarity"
                                write (output_unit, *) " "//reci//"; reenter ..."
                                exit
                            end if
                            nkp = nkp + 1
                            if (nkp > ndkp) then
                                call ralloc (kp, nkp, "KP", "GETRMP" )
                                ndkp = nkp
                            end if
                            if (lelec) then
                                kp(nkp) = mult*(-1)**mult
                            else if (lmagn) then
                                kp(nkp) = mult*(-1)**(mult + 1)
                            end if
                            istart = 0
                        end if
                    end if
                end do
                ! An invalid token stops the scan and retries the input.
                if (i > len(record)) exit
            end do
            if (nkp /= 0) exit
        end do
    end if !myid=0
    call mpi_bcast(nkp,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    !
    !   Trim array to the exact size
    !
    if (ndkp /= nkp) call ralloc (kp, nkp, "KP", "GETRMP")
    call mpi_bcast(kp(1),nkp,mpi_integer,0,mpi_comm_world, mpi_error_code)
    !
    !   If M1 or E2 inquire if the transitions are between levels
    !   with different J quantum numbers.
    !
    if(mpi_rank == 0) then
        do i = 1, nkp
            if (kp(i) == 1) then
                write (*, *)                                         &
                    "M1 transitions only between levels with different J?"
                yes = getyn()
                if (yes) then
                    noffd1 = 1
                else
                    noffd1 = 0
                end if
            end if
            if (kp(i) /= 2) cycle
            write (*, *)                                            &
                "E2 transitions only between levels with different J?"
            yes = getyn()
            if (yes) then
                noffd2 = 1
            else
                noffd2 = 0
            end if
        end do
    end if !myid=0
    call mpi_bcast(noffd1,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
    call mpi_bcast(noffd2,1,mpi_integer,0,mpi_comm_world, mpi_error_code)
    !
    return
end subroutine read_multipole_options
