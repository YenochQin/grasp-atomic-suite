!***********************************************************************
!                                                                      *
subroutine configure_debug_output
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    use debug_output_settings, only: angular_debug_unit
    !-----------------------------------------------
    !                                                                      *
    !   This subroutine sets the arrays that control debug printout from   *
    !   the radial and angular modules of the GRASP92 suite.               *
    !                                                                      *
    !   Call(s) to: [LIB92]: GETYN, LENGTH, OPENFL.                        *
    !                                                                      *
    !   Written by Farid A Parpia               Last update: 24 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    !...Translated by Gediminas Gaigalas 11/18/19
    use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
    use calculation_input_defaults, only: inpci, ndef, ndump
    use yes_no_input, only: getyn
    use file_open, only: openfl
    implicit none (type, external)

    ! Local declarations
    integer   :: open_status
    logical   :: yes
    character :: file_path*256, default_file_path*11, form*11, status*3
    !-----------------------------------------------
    !
    !   Initialise the arrays that control the debug printout
    !
    ldbpa = .false.
    !
    ldbpg = .false.
    !
    ldbpr = .false.
    !
    if (ndef == 0) return

    write (output_unit, *) "Generate debug printout?"
    yes = getyn()
    if (yes) then
        !
        !   The  .dbg  file is formatted; open it on unit 99
        !
        default_file_path = "rdensity.dbg"
        form = "FORMATTED"
        status = "NEW"
        !
        write (output_unit, *) "File  rdensity.dbg  will be created as the"
        write (output_unit, *) " rdensity DeBuG Printout File; enter another"
        write (output_unit, *) " file name if this is not acceptable;"
        write (output_unit, *) " null otherwise:"
        read (*, "(A)") file_path
        !
        if (len_trim(file_path) == 0) file_path = default_file_path
        !
        do
            call openfl (angular_debug_unit, file_path, form, status, open_status)
            if (open_status == 0) exit
            do
                write (output_unit, *) "Enter a name for the rdensity DeBuG Printout"
                write (output_unit, *) " file that is to be created:"
                read (*, "(A)") file_path
                if (len_trim(file_path) /= 0) exit
            end do
        end do
        !
        !   Set options for general printout
        !
        write (output_unit, *) " Print out the machine constants used?"
        yes = getyn()
        if (yes) ldbpg(1) = .true.
        write (output_unit, *) " Print out the physical constants used?"
        yes = getyn()
        if (yes) ldbpg(2) = .true.
        !
        !   Set options for radial modules
        !
        write (output_unit, *) " Printout from radial modules?"
        yes = getyn()
        if (yes) then
            write (output_unit, *) " Printout from RADGRD?"
            yes = getyn()
            if (yes) ldbpr(1) = .true.
            write (output_unit, *) " Printout from NUCPOT?"
            yes = getyn()
            if (yes) ldbpr(2) = .true.
            write (output_unit, *) " Printout from LODRWF?"
            yes = getyn()
            if (yes) ldbpr(3) = .true.
            !
        end if
        !
        !   Set options for angular modules
        !
        write (output_unit, *) " Printout from angular modules?"
        yes = getyn()
        if (yes) then
            write (output_unit, *) " Printout from LODCSL?"
            yes = getyn()
            if (yes) ldbpa(1) = .true.
            write (output_unit, *) " Print out V coefficients?"
            yes = getyn()
            if (yes) ldbpa(2) = .true.
        end if
        !
    end if
    !
    return
end subroutine configure_debug_output
