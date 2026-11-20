!***********************************************************************
!                                                                      *
!     This program transforms the initial and final state radial       *
!     orbitals to a biorthonormal set. The CI coefficients are         *
!     then counter transformed as to leave both the initial and        *
!     final state wave functions invariant                             *
!                                                                      *
!     General references:                                              *
!                                                                      *
!     P.A. Malmqvist, Int.J. of Quantum Chemistry, XXX, 479-94 (1986)  *
!     J. Olsen, M.R. Godefroid, P. Jonsson, P.A. Malmqvist and         *
!     C. Froese Fischer, Phys. Rev. E, 4499 (1995)                     *
!                                                                      *
program rbiotransform_main
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: output_unit
    !                                                                      *
    !     Program written by                                               *
    !                                                                      *
    !     Per Jonsson, Department of Computer Science                      *
    !                  Vanderbilt University, USA                          *
    !                                                                      *
    !     Date June 1996                                                   *
    !                                                                      *
    !     Modified by Gediminas Gaigalas for new spin-angular integration  *
    !     and for reducing usage of CPU memory.        NIST, October 2017  *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    use grasp_parameters,   only: max_orbitals, max_radial_points
    use calculation_input_defaults,       only: ndef, ndump
    use biorthogonal_subshell_layout,         only: nlmax, kamax, nshlii, nshlff
    use biorthogonal_orbital_sets,         only: nwii, nwff
    use radial_wavefunctions,          only: pfii, qfii, pfff, qfff
    use biorthogonal_orbital_rotations,         only: cici, cfci
    use yes_no_input, only: getyn
    use machine_constants, only: setmc
    use physical_constants, only: setcon
    use isotope_data_input, only: setiso

    ! Procedure interfaces
    use configure_radial_grid_interface, only: configure_radial_grid

    ! Module dependencies
    use radial_grid_setup, only: radgrd
    use radial_stencil_weights, only: setqic

    ! Procedure interfaces
    use read_state_filenames_interface, only: read_state_filenames
    use load_merged_csf_basis_interface, only: load_merged_csf_basis
    use copy_state_csf_data_interface, only: copy_state_csf_data
    use prepare_biorthogonal_subshells_interface, only: prepare_biorthogonal_subshells
    use read_initial_radial_wavefunctions_interface, only: read_initial_radial_wavefunctions
    use read_final_radial_wavefunctions_interface, only: read_final_radial_wavefunctions
    use calculate_orbital_overlap_matrix_interface, only: calculate_orbital_overlap_matrix
    use calculate_subshell_overlap_interface, only: calculate_subshell_overlap
    use transform_orbital_sets_interface, only: transform_orbital_sets
    use write_transformed_radial_wavefunctions_interface, only: write_transformed_radial_wavefunctions
    use prepare_biorthogonal_mcp_interface, only: prepare_biorthogonal_mcp

    ! Module dependencies
    use run_timer_start, only: starttime
    use run_timer_finish, only: stoptime
    implicit none (type, external)
    !GG      INTEGER, PARAMETER :: LWORK1 = 100000
    integer, parameter :: lwork1 = 10000000
    integer :: ninshli(nlmax), ninshlf(nlmax)
    integer :: ntestg, ntestl, ntest, inpci, ncore1, ncore2, mxl
    integer :: start_clock_count
    real(real64) :: work(lwork1)
    real(real64) :: cishl(nlmax*nlmax)
    real(real64)   :: s(max_orbitals, max_orbitals)
    real(real64) :: cfshl(nlmax*nlmax)
    logical :: yes
    character :: name(2)*24
    character(len=128) :: isotope_file_path
    !-----------------------------------------------
    !
    isotope_file_path = "isodata"
    !   Debug flags
    !
    ntestg = 0
    ntestl = 0
    ntest = max0(ntestl,ntestg)
    !
    call starttime (start_clock_count, "BIOTRANSFORM")
    write(*,*)
    write(*,*) "RBIOTRANSFORM"
    write(*,*) "This program transforms the initial and final wave"
    write(*,*) "functions so that standard tensor albegra can be"
    write(*,*) "used in evaluation of the transition parameters"
    write(*,*) "Input files:  isodata, name1_descriptors.parquet, name1.wfn, name1.mix/name1.cimix"
    write(*,*) "              name2_descriptors.parquet, name2.wfn, name2.mix/name2.cimix         "
    print*,"              name1.TB, name2.TB (optional angular files)"
    write(*,*) "Output files: name1.bw, name1.(c)bm, "
    write(*,*) "              name2.bw, name2.c(bm)"
    write(*,*) "              name1.TB, name2.TB (angular files)"
    write(*,*)

    write (output_unit, *) "Default settings?"
    yes = getyn()
    write (output_unit, *)
    if (yes) then
        ndef = 0
        ndump = 1
    else
        ndef = 1
        write (output_unit, *) "Dump angular data on file?"
        yes = getyn()
        write (output_unit, *)
        if (yes) then
            ndump = 1
        else
            ndump = 0
        end if
    end if
    write (output_unit, *) "Input from a CI calculation?"
    yes = getyn()
    write (output_unit, *)
    if (yes) then
        inpci = 0
    else
        inpci = 1
    end if
    !
    !   Perform machine- and installation-dependent setup
    !
    call setmc
    !
    !   Set up the physical constants
    !
    call setcon
    !
    !   Open, check, load data from, and close the  .iso  file
    !
    call setiso (isotope_file_path)
    !
    !   Determine the parameters controlling the radial grid
    !
    call configure_radial_grid
    !
    !   Generate the radial grid
    !
    call radgrd
    !
    !   Set up the coefficients for the numerical procedures
    !
    call setqic
    !
    !   Obtain the names of the initial and final state files
    !   and open files where the transformed orbitals and CI
    !   coefficients are to be dumped
    !
    call read_state_filenames (name)
    !
    !   Open, check, load data from and close the initial state CSL file.
    !
    call load_merged_csf_basis (name(1), ncore1,1)
    !
    !   Transfer the data to the initial state COMMON
    !
    call copy_state_csf_data (1)
    !
    !   Open, check, load data from and close the final state CSL file.
    !
    call load_merged_csf_basis (name(2), ncore2,2)
    !
    !   Transfer the data to the final state COMMON
    !
    call copy_state_csf_data (2)
    !
    !   Determine the number of kappa quantum numbers and
    !   the number of orbitals for each kappa quantum number
    !   for the initial state and final states
    !
    call prepare_biorthogonal_subshells (ntestg, ncore1, ncore2)
    !
    !   Read the the radial orbitals for the initial state
    !
    call read_initial_radial_wavefunctions (name(1), ntestg)
    !
    !   Read the the radial orbitals for the initial state
    !
    call read_final_radial_wavefunctions (name(2), ntestg)
    !
    !   Calculate the radial overlap matrices
    !
    write (*, *)
    write (*, *) " ******************************************"
    write (*, *) "  Overlap matrix before orbital rotations"
    write (*, *) " *****************************************"
    write (*, *)

    call calculate_orbital_overlap_matrix

    call calculate_subshell_overlap (s, nwii, nwff)

    !
    !   Once we have the overlap matrices
    !   we can manipulate the initial and final state separately.
    !
    mxl = kamax
    !
    !. Calculate biorthonormal orbitals, and orbital matrix
    !. for counter transformation of CI coefficients.
    !
    call transform_orbital_sets (pfii, qfii, nshlii, ninshli, pfff, qfff, &
        nshlff, ninshlf, max_radial_points, kamax, work, lwork1, ntestg, &
        cishl, cici, cfshl, cfci)
    write (*, *)
    write (*, *) " ****************************************"
    write (*, *) "  Overlap matrix after orbital rotations"
    write (*, *) " ****************************************"
    write (*, *)

    call calculate_orbital_overlap_matrix
    !
    !  Write the transformed radial functions to file
    !
    call write_transformed_radial_wavefunctions (name)
    !
    !   Obtain one-electron coupling coefficients for the initial state.
    !   The coefficients are dumped on files one kappa in turn and
    !   thus the different kappa can be manipulated independently.
    !   The interface with the transformation part is in the routine mcp
    !
    call prepare_biorthogonal_mcp (name(1), 1, ntestg, inpci)
    !
    !   Obtain one-electron coupling coefficients for the final state.
    !   The coefficients are dumped on files one kappa in turn and
    !   thus the different kappa can be manipulated independently.
    !   The interface with the transformation part is in the routine mcp
    !
    call prepare_biorthogonal_mcp (name(2), 2, ntestg, inpci)

    call stoptime (start_clock_count, "BIOTRANSFORM")

    stop
end program rbiotransform_main
