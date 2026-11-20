!***********************************************************************
!                                                                      *
subroutine connect_ci_restart_file(isotope_file_path, radial_file_path, block_labels)
    ! Module dependencies
    !                                                                      *
    !   Open, check, load data from the  .res  file.                       *
    !                                                                      *
    !   Call(s) to: [LIB92]: GETYN, OPENFL.
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
    !   Modified by Xinghong                  Last revision: 23 Jun 1998   *
    !                                                                      *
    !***********************************************************************
    use array_memory_management, only: alloc
    use calculation_input_defaults, only: inpci, ndef, ndump
    use ci_restart_units, only: restart_matrix_unit
    use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
    use csf_block_cutoffs, only: iccut, iccutblk
    use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
    use yes_no_input, only: getyn
    use file_open, only: openfl

    ! Procedure interfaces
    use read_ci_restart_interface, only: read_ci_restart
    use read_ci_options_interface, only: read_ci_options
    implicit none (type, external)

    ! Arguments
    character  :: isotope_file_path*(*)
    character  :: radial_file_path*(*)
    character(len=8)  :: block_labels(*)

    ! Local constants
    character(len=11), parameter :: form = "UNFORMATTED"
    character(len=7), parameter :: status = "UNKNOWN"
    character(len=6), parameter :: restitle = "R92RES"

    ! Local declarations
    integer :: open_status, io_status
    logical :: found, is_restart
    character :: r92res*6, file_path*11, block_label_text*3
    !-----------------------------------------------
    !
    ! Compose the "rci.res" file name
    !
    file_path = "rci.res"
    !
    ! Ask if this is a restart
    !
    if (ndef /= 0) then
        write (stderr_unit, *) "Restarting RCI90 ?"
        is_restart = getyn()
    else
        is_restart = .false.
    end if
    !      IF (RESTRT) THEN
    !        WRITE(734,'(a)') 'y            ! Restarting RCI90 ?'
    !      ELSE
    !        WRITE(734,'(a)') 'n            ! Restarting RCI90 ?'
    !      END IF
    !
    ! Do some settings and checks
    !
    if (is_restart) then
        !         ...Restart, make sure file exist
        inquire(file=file_path, exist=found)
        if (.not.found) stop "setres: .res does not exist"
    end if
    !
    ! Open the .res file
    !
    call openfl (restart_matrix_unit, file_path, form, status, open_status)
    if (open_status /= 0) stop "setres: Error openning .res file"
    !
    ! If restart, load the contents. Otherwise generate them via getcid
    !
    ! But first of all, iccutblk() is needed in both cases
    !
    call alloc (iccutblk, symmetry_block_count, "ICCUTBLK", "SETRES")

    if (is_restart) then
        !        ...Check the signature of the file
        read (restart_matrix_unit, iostat=io_status) r92res
        if (io_status/=0 .or. r92res/=restitle) then
            close(restart_matrix_unit)
            stop "setres: Not RCI92 .res file"
        end if

        !         ...Read and check restart information
        call read_ci_restart

    else

        !         ...Write the file header
        !         ...Generate the first part of the .res file
        write (restart_matrix_unit) restitle
        call read_ci_options (isotope_file_path, radial_file_path, block_labels)

    end if

    return
end subroutine connect_ci_restart_file
