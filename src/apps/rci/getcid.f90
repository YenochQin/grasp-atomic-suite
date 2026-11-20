!***********************************************************************
!                                                                      *
      subroutine getcid (isotope_file_path, radial_file_path, block_labels)
      use calculation_file_units, only: calculation_log_unit
!                                                                      *
!   Interactively determines the data governing the CI problem.        *
!   iccut is replaced by an array iccutblk(1:nblock)
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, NUCPOT, RADGRD, SETQIC.                *
!               [RCI92]: SETISO, SETRWF.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use atomic_calculation_state,           only: &
          vacuum_speed_of_light_au, speed_of_light_au, nuclear_charge, accy, electron_count, emn
      use calculation_input_defaults, only: inpci, ndef, ndump
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use csf_block_cutoffs, only: iccut, iccutblk
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_potential_state,          only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
      use isotope_properties, only: dmomnm, hfsd, hfsi, hfsq, qmomb, smsd, smsi, smsq, sqn
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_energy_parameters, &
          orbital_origin_exponents
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use breit_frequency_scale, only: wfact
      use rci_prerun_settings, only: coeffcut1, coeffcut2, iprerun, ncsfpre
      use where_c, only: restart_matrix_unit
      use qed_orbital_cutoffs, only: nqedcut, nqedmax
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use nuclear_coulomb_potential, only: nucpot
      use radial_grid_setup, only: radgrd
      use isotope_data_input, only: setiso
      use radial_stencil_weights, only: setqic
      use radial_wavefunction_input, only: setrwfa
      implicit none (type, external)
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      character(len=*) :: isotope_file_path, radial_file_path
      character(len=8) ::  block_labels(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       logical :: yes
       integer :: symmetry_block_index, jblock, ntmp, np10, orbital_index
!-----------------------------------------------------------------------
!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole
! common block.
!
!      PRINT *, 'Calling SETISO ...'
      call setiso (isotope_file_path)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
      if (ndef /= 0) then
         write (stderr_unit,*) "Revise the physical speed of light (",vacuum_speed_of_light_au,   &
                           " in a.u.) ?"
         yes = getyn ()
         if (yes) then
            write (stderr_unit,*) "Enter the revised value:"
            read *,speed_of_light_au
         else
            speed_of_light_au = vacuum_speed_of_light_au
         end if
      else
         speed_of_light_au = vacuum_speed_of_light_au
      end if
!
! Treat some CSF's as perturbation ? Broadcasting is used only in
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
      if (ndef /= 0) then
         write (stderr_unit,*) "Treat contributions of some CSFs",            &
                           " as first-order perturbations?"
         lfordr = getyn ()
      else
         lfordr = .false.
      end if

! Get iccutblk() from the user-input

      if (.not. lfordr) then
         !...Default first
         do symmetry_block_index = 1, symmetry_block_count
           iccutblk(symmetry_block_index) = csf_counts_by_block(symmetry_block_index)
         end do
      else

      ! Let master do the i/o, then broadcast
            write (stderr_unit,*) "There are ", symmetry_block_count, "blocks. They are:"
            write (stderr_unit,*) "  block     J Parity     No of CSFs"
            do symmetry_block_index = 1, symmetry_block_count
               write (stderr_unit,*) symmetry_block_index, block_labels(symmetry_block_index)(1:5), &
                  & csf_counts_by_block(symmetry_block_index)
            end do

            write (stderr_unit,*)
            write (stderr_unit,*) "Enter iccut for each block"
            do jblock = 1, symmetry_block_count
               write (stderr_unit,*) "Block ", jblock, "   ncf = ",           &
                              csf_counts_by_block(jblock)                            &
                      , " id = ", block_labels(jblock)(1:5)
  123          read (stdin_unit,*) ntmp
               if (ntmp >= 0 .and. ntmp <= csf_counts_by_block(jblock)) then
                  iccutblk(jblock) = ntmp
               else
                  write (stderr_unit,*) "ICCUT out of range, re-enter:"
                  goto 123
               end if
               write(calculation_log_unit,*) ntmp,"! ICCUT for block",jblock
            end do
      end if

!*****************************************************************
!
! Pre-run ?
!
!     IF (IPRERUN .EQ. 0) THEN

!        WRITE (istde,*) ' Prerun with limited interaction?'
!        YES = GETYN ()

!        IF (YES) THEN
!           IPRERUN = 1
!           LTRANS = .FALSE.
!           LVP = .FALSE.
!           LNMS = .FALSE.
!           LSMS = .FALSE.
!           LSE = .FALSE.

!           WRITE (istde,*)  ' Give CSL cut'
!           READ *, NCSFPRE
!           WRITE (istde,*)  ' Give coefficient cut for H_0'
!           READ *, COEFFCUT1
!           WRITE (istde,*) ' Give coefficient cut for the transvers'
!    &,                  ' interaction'
!           GOTO 99
!        ENDIF
!     ENDIF
!*****************************************************************
!
! Include transverse ?
! Quantities to be obtained: LTRANS, WFACT
!
      write (stderr_unit,*) "Include contribution of H (Transverse)?"
      ltrans = getyn ()
      write (stderr_unit,*) "Modify all transverse photon frequencies?"
      yes = getyn ()
      if (yes) then
         write (stderr_unit,*) "Enter the scale factor:"
         read *, wfact
      else
         wfact = 1.0_real64
      end if
!
! Other interactions ? One logical for each case. Done altogether
!
      write (stderr_unit,*) "Include H (Vacuum Polarisation)?"
      lvp = getyn ()

      write (stderr_unit,*) "Include H (Normal Mass Shift)?"
      lnms = getyn ()

      write (stderr_unit,*) "Include H (Specific Mass Shift)?"
      lsms = getyn ()

      write (stderr_unit,*) "Estimate self-energy?"
      lse = getyn ()
      if (lse.eqv..true.) then
         nqedcut = 1
         write (stderr_unit,*)                                             &
      "Largest n quantum number for including self-energy for orbital"
         write (stderr_unit,*) "n should be less or equal 8"
         read *, nqedmax
      else
         nqedcut = 0
      end if

      if (ltrans) then
        write(calculation_log_unit,"(a)") "y            ! Contribution of H Transverse?"
      else
        write(calculation_log_unit,"(a)") "n            ! Contribution of H Transverse?"
      end if
      if (yes) then
        write(calculation_log_unit,"(a)") "y            ! Modify photon frequencies?"
        write(calculation_log_unit,*) wfact,"! Scale factor"
      else
        write(calculation_log_unit,"(a)") "n            ! Modify photon frequencies?"
      end if
      if (lvp) then
        write(calculation_log_unit,"(a)") "y            ! Vacuum polarization?"
      else
        write(calculation_log_unit,"(a)") "n            ! Vacuum polarization?"
      end if
      if (lnms) then
        write(calculation_log_unit,"(a)") "y            ! Normal mass shift?"
      else
        write(calculation_log_unit,"(a)") "n            ! Normal mass shift?"
      end if
      if (lsms) then
        write(calculation_log_unit,"(a)") "y            ! Specific mass shift?"
      else
        write(calculation_log_unit,"(a)") "n            ! Specific mass shift?"
      end if
      if (lse) then
        write(calculation_log_unit,"(a)") "y            ! Self energy?"
        write(calculation_log_unit,*) nqedmax, "! Max n for including self energy"
      else
        write(calculation_log_unit,"(a)") "n            ! Self energy?"
      end if
!
! Parameters controlling the radial grid
!
!cff-cjb grid
!cff use default grid
!
! 99  IF (NPARM .EQ. 0) THEN
!        RNT = EXP (-65.D0/16.D0) / Z
!        H = 0.5D0**4
!        N = MIN (220,NNNP)
!     ELSE
!CFF     .. should be Z-dependent
         radial_grid_scale = 2.0e-6_real64/nuclear_charge
         transformed_grid_step = 5.0e-2_real64
         radial_point_count = max_radial_points
!     ENDIF
      asymptotic_grid_step = 0.0_real64
!     IF ( NDEF.NE.0) THEN
!        WRITE (istde,*) 'The default radial grid parameters',         &
!                       ' for this case are:'
!        WRITE (istde,*) ' RNT = ',RNT,';'
!        WRITE (istde,*) ' H = ',H,';'
!        WRITE (istde,*) ' HP = ',HP,';'
!        WRITE (istde,*) ' N = ',N,';'
!        WRITE (istde,*) ' revise these values?'
!        YES = GETYN ()
!
!        IF (YES) THEN
!           WRITE (istde,*) 'Enter RNT:'
!           READ *, RNT
!           WRITE (istde,*) 'Enter H:'
!           READ *, H
!           WRITE (istde,*) 'Enter HP:'
!           READ *, HP
!           WRITE (istde,*) 'Enter N:'
!           READ *, N
!        ENDIF
!     ENDIF
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = transformed_grid_step**6
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Generate $- r \times V_ (r)$
!
      call nucpot
!
!   Load the radial wavefunctions
!
      call setrwfa (radial_file_path)
!
!   Write the basic parameters of the model electron cloud to the
!   .res  file; this is the second record on the file --- the
!   first is the header (see SUBROUTINE SETRES)
!
      write (restart_matrix_unit) electron_count, csf_count, orbital_count, symmetry_block_count ! ncf is ncftot, from setcsll
!
!   Write the nuclear parameters and $- r \times V_ (r)$
!   to the  .res  file; these are the third, fourth, and fifth
!   records on the file
!
      write (restart_matrix_unit) nuclear_charge, emn
      write (restart_matrix_unit) nuclear_model_parameter_count,(nuclear_model_parameters(symmetry_block_index), &
         & symmetry_block_index = 1, nuclear_model_parameter_count)
      write (restart_matrix_unit) radial_point_count, (negative_radius_times_potential(symmetry_block_index), &
         & symmetry_block_index = 1, radial_point_count), &
         & first_grid_index_beyond_fermi_c
!
!   Write the physical effects specification to the  .res  file.
!   iccutblk() is now an array of length nblock.
!   This is the sixth record on the file
!
      write (restart_matrix_unit) speed_of_light_au, lfordr, (iccutblk(symmetry_block_index), symmetry_block_index = 1, &
         & symmetry_block_count),            &
                    ltrans, wfact, lvp, lnms, lsms
!
!   Write the grid data to the  .res  file; this is the seventh
!   record on the file
!
      np10 = radial_point_count + 10
      write (restart_matrix_unit) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, &
         & (radial_grid(symmetry_block_index), symmetry_block_index = 1, np10), &
         & &
                 (radial_grid_derivative(symmetry_block_index), symmetry_block_index = 1, np10), &
                    & (radial_derivative_over_radius(symmetry_block_index), symmetry_block_index = 1, np10) ! (restart_matrix_unit)
!
!   Write out the interpolated radial wavefunctions; there are
!   2*NW such records; thus the total number of records written
!   at this point is 7+2*NW
!
      do orbital_index = 1, orbital_count
         write (restart_matrix_unit) orbital_energy_parameters(orbital_index), orbital_origin_exponents(orbital_index), &
            & orbital_origin_parameters(orbital_index), &
            & orbital_last_grid_index(orbital_index)
         write (restart_matrix_unit) (radial_large_components(symmetry_block_index,orbital_index), symmetry_block_index = 1, &
            & orbital_last_grid_index(orbital_index)), &
            (radial_small_components(symmetry_block_index,orbital_index) &
            & , symmetry_block_index = 1, orbital_last_grid_index(orbital_index))
      end do
!
      return
      end subroutine getcid
