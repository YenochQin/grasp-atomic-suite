

!***********************************************************************
!***********************************************************************
!***********************************************************************
!**                                                                  ***
!**               *******  *******   **    **  *******               ***
!**               **       **    **  **    **  **                    ***
!**               **       **    **  **    **  **                    ***
!**               *****    *******   ** ** **  *****                 ***
!**               **       **  **    ** ** **  **                    ***
!**               **       **   **   ** ** **  **                    ***
!**               *******  **    **   **  **   **                    ***
!**                                                                  ***
!**   Program to create a  .rwf  file for RSCF92 by merging one or   ***
!**   more existing  .rwf  files or generating Thomas-Fermi or hy-   ***
!**   drogenic functions.                                            ***
!**                                                                  ***
!**                              GRASP92                             ***
!**           F. A. Parpia, C. F. Fischer, and I. P. Grant           ***
!**                                                                  ***
!**                                                                  ***
!***********************************************************************
!***********************************************************************
!***********************************************************************
!                                                                      *
program rwfnestimate_main
    ! Module dependencies
    !                                                                      *
    !   Entry routine for RCI92. Controls the entire computation.          *
    !                                                                      *
    !   Call(s) to: [LIB92]: SETMC, SETCON.                                *
    !               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
    !                        SETRES, SETSUM, STRSUM.                       *
    !               [NJGRAF]: FACTT.                                       *
    !                                                                      *
    !   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
    !                                                                      *
    !***********************************************************************
    !
    use calculation_input_defaults, only: nondefault_options_flag
    use standard_io_units, only: stderr_unit
    use rwfnestimate_input, only: radial_input_options, parse_radial_options, load_radial_header
    use calculation_state_files, only: check_orbital_output
    use calculation_state_files, only: state_files, state_name_length, select_state, files_for_state, orbital_input_path

    ! Procedure interfaces
    use configure_debug_output_interface, only: configure_debug_output

    ! Module dependencies
    use machine_constants, only: initialize_machine_constants
    use physical_constants, only: initialize_physical_constants

    ! Procedure interfaces
    use open_summary_file_interface, only: open_summary_file
    use set_hydrogenic_screening_interface, only: set_hydrogenic_screening
    use read_calculation_options_interface, only: read_calculation_options
    use write_input_summary_interface, only: write_input_summary
    use generate_initial_radial_wavefunctions_interface, only: generate_initial_radial_wavefunctions

    ! Module dependencies
    use radial_orbital_orthogonalization, only: orthogonalize_radial_orbitals

    ! Procedure interfaces
    use write_initial_radial_wavefunctions_interface, only: write_initial_radial_wavefunctions
    implicit none (type, external)
    type(state_files) :: files
    character(state_name_length) :: state
    integer :: core_subshell_count
    logical :: yes
    type(radial_input_options) :: input_options
    !-----------------------------------------------

    !   Startup message
    !
    call parse_radial_options(input_options)
    write (stderr_unit, *) "RWFNESTIMATE"
    write (stderr_unit, *) "This program estimates radial wave functions"
    write (stderr_unit, *) "for orbitals"
    write (stderr_unit, *) "Input files: isodata, descriptor header TOML, optional rwfn file"
    call select_state(state, "Default settings ?", yes)
    files = files_for_state(state)
    call check_orbital_output(files%initial_orbitals, writing_initial=.true.)
    if (input_options%filename == "") input_options%filename = files%header
    block
        ! Local declarations
        character(:), allocatable :: basename
        basename = input_options%filename(scan(input_options%filename, "/"//achar(92), back=.true.)+1:)
        if (basename == files%initial_orbitals .or. basename == files%estimate_debug .or. &
            basename == files%estimate_summary) error stop "Header uses a reserved output name"
    end block
    write (stderr_unit, *) "Selected input: ", input_options%filename
    write (stderr_unit, *) "Output file: ", files%initial_orbitals


    !
    if (yes) then
        nondefault_options_flag = 0
    else
        nondefault_options_flag = 1
    end if

    !
    !   Determine if there is to be any debug printout; this will be
    !   made on the  .dbg  file
    !
    call configure_debug_output(files%estimate_debug)
    !
    !   Perform machine- and installation-dependent setup
    !
    call initialize_machine_constants
    !
    !   Set up the physical constants
    !
    call initialize_physical_constants
    !
    !   Open the  .sum  file
    !
    if (nondefault_options_flag /= 0) call open_summary_file(files%estimate_summary)
    !
    !   Load orbital metadata from the selected rCSFs header.
    !
    call load_radial_header(input_options, core_subshell_count)
    !
    !   Hydrogenic screen parameters for all orbitals
    !
    call set_hydrogenic_screening (core_subshell_count)
    !
    !   Determine other relevant information
    !
    call read_calculation_options
    !
    !   Write the first part of the  .sum  file
    !
    if (nondefault_options_flag /= 0) call write_input_summary
    !
    !   Generate the subshell radial wavefunctions
    !
    if (state /= "") files%orbitals = orbital_input_path(state)
    call generate_initial_radial_wavefunctions(files%orbitals)
    !
    !   Orthogonalize the radial orbitals
    !
    call orthogonalize_radial_orbitals
    !
    !   Write the subshell radial wavefunctions out
    !
    call write_initial_radial_wavefunctions(files%initial_orbitals)
    !
    !   Print completion message
    !
    write (stderr_unit, *) "RWFNESTIMATE: Execution complete."
    !
    stop
end program rwfnestimate_main
