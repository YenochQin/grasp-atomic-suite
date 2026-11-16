!***********************************************************************
!                                                                      *
      subroutine getcid (isofile, rwffile, idblk)
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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_radial_points
      use decide_c
      use def_c,           only: vacuum_speed_of_light_au, speed_of_light_au, nuclear_charge, accy, electron_count, emn
      use default_c
      use grid_c
      use hblock_c
      use iccu_c
      use iounit_c
      use npar_c
      use npot_c,          only: negative_radius_times_potential, first_grid_index_beyond_fermi_c
      use nsmdat_c
      use orb_c
      use wave_c
      use wfac_c
      use blim_c
      use where_c
      use qedcut_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use nucpot_i, only: nucpot
      use radgrd_i, only: radgrd
      use setiso_i, only: setiso
      use setqic_i, only: setqic
      use setrwfa_i, only: setrwfa
      implicit none (type, external)
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      character(len=*) :: isofile, rwffile
      character(len=8) ::  idblk(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       logical :: yes
       integer :: i, jblock, ntmp, np10, j
!-----------------------------------------------------------------------
!
! Open, check, load data from, and close the  .iso  file
! Data loaded are: EMN,Z,/NPAR/,/NSMDAT/ where /.../ means whole
! common block.
!
!      PRINT *, 'Calling SETISO ...'
      call setiso (isofile)
!
! The speed of light, if non-default then spread from node-0
! Quantities to be obtained: C
!
      if (ndef .ne. 0) then
         write (stderr_unit,*) 'Revise the physical speed of light (',vacuum_speed_of_light_au,   &
                           ' in a.u.) ?'
         yes = getyn ()
         if (yes) then
            write (stderr_unit,*) 'Enter the revised value:'
            read *,speed_of_light_au
         else
            speed_of_light_au = vacuum_speed_of_light_au
         endif
      else
         speed_of_light_au = vacuum_speed_of_light_au
      endif
!
! Treat some CSF's as perturbation ? Broadcasting is used only in
! non-default mode, as the above case for speed of light.
! Quantities to be obtained: LFORDR
!
      if (ndef .ne. 0) then
         write (stderr_unit,*) 'Treat contributions of some CSFs',            &
                           ' as first-order perturbations?'
         lfordr = getyn ()
      else
         lfordr = .false.
      endif

! Get iccutblk() from the user-input

      if (.not. lfordr) then
         !...Default first
         do i = 1, symmetry_block_count
           iccutblk(i) = csf_counts_by_block(i)
         enddo
      else

      ! Let master do the i/o, then broadcast
            write (stderr_unit,*) 'There are ', symmetry_block_count, 'blocks. They are:'
            write (stderr_unit,*) '  block     J Parity     No of CSFs'
            do i = 1, symmetry_block_count
               write (stderr_unit,*) i, idblk(i)(1:5), csf_counts_by_block(i)
            enddo

            write (stderr_unit,*)
            write (stderr_unit,*) 'Enter iccut for each block'
            do jblock = 1, symmetry_block_count
               write (stderr_unit,*) 'Block ', jblock, '   ncf = ',           &
                              csf_counts_by_block(jblock)                            &
                      , ' id = ', idblk(jblock)(1:5)
  123          read (stdin_unit,*) ntmp
               if (ntmp .ge. 0 .and. ntmp .le. csf_counts_by_block(jblock)) then
                  iccutblk(jblock) = ntmp
               else
                  write (stderr_unit,*) 'ICCUT out of range, re-enter:'
                  goto 123
               endif
               write(734,*) ntmp,'! ICCUT for block',jblock
            enddo
      endif

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
      write (stderr_unit,*) 'Include contribution of H (Transverse)?'
      ltrans = getyn ()
      write (stderr_unit,*) 'Modify all transverse photon frequencies?'
      yes = getyn ()
      if (yes) then
         write (stderr_unit,*) 'Enter the scale factor:'
         read *, wfact
      else
         wfact = 1.0_real64
      endif
!
! Other interactions ? One logical for each case. Done altogether
!
      write (stderr_unit,*) 'Include H (Vacuum Polarisation)?'
      lvp = getyn ()

      write (stderr_unit,*) 'Include H (Normal Mass Shift)?'
      lnms = getyn ()

      write (stderr_unit,*) 'Include H (Specific Mass Shift)?'
      lsms = getyn ()

      write (stderr_unit,*) 'Estimate self-energy?'
      lse = getyn ()
      if (lse.eqv..true.) then
         nqedcut = 1
         write (stderr_unit,*)                                             &
      'Largest n quantum number for including self-energy for orbital'
         write (stderr_unit,*) 'n should be less or equal 8'
         read *, nqedmax
      else
         nqedcut = 0
      end if

      if (ltrans) then
        write(734,'(a)') 'y            ! Contribution of H Transverse?'
      else
        write(734,'(a)') 'n            ! Contribution of H Transverse?'
      end if
      if (yes) then
        write(734,'(a)') 'y            ! Modify photon frequencies?'
        write(734,*) wfact,'! Scale factor'
      else
        write(734,'(a)') 'n            ! Modify photon frequencies?'
      end if
      if (lvp) then
        write(734,'(a)') 'y            ! Vacuum polarization?'
      else
        write(734,'(a)') 'n            ! Vacuum polarization?'
      end if
      if (lnms) then
        write(734,'(a)') 'y            ! Normal mass shift?'
      else
        write(734,'(a)') 'n            ! Normal mass shift?'
      end if
      if (lsms) then
        write(734,'(a)') 'y            ! Specific mass shift?'
      else
        write(734,'(a)') 'n            ! Specific mass shift?'
      end if
      if (lse) then
        write(734,'(a)') 'y            ! Self energy?'
        write(734,*) nqedmax, '! Max n for including self energy'
      else
        write(734,'(a)') 'n            ! Self energy?'
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
         transformed_grid_step = 5.e-2_real64
         radial_point_count = max_radial_points
!     ENDIF
      asymptotic_grid_step = 0._real64
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
      call setrwfa (rwffile)
!
!   Write the basic parameters of the model electron cloud to the
!   .res  file; this is the second record on the file --- the
!   first is the header (see SUBROUTINE SETRES)
!
      write (imcdf) electron_count, csf_count, orbital_count, symmetry_block_count ! ncf is ncftot, from setcsll
!
!   Write the nuclear parameters and $- r \times V_ (r)$
!   to the  .res  file; these are the third, fourth, and fifth
!   records on the file
!
      write (imcdf) nuclear_charge, emn
      write (imcdf) nuclear_model_parameter_count,(nuclear_model_parameters(i), i = 1, nuclear_model_parameter_count)
      write (imcdf) radial_point_count, (negative_radius_times_potential(i), i = 1, radial_point_count), &
         & first_grid_index_beyond_fermi_c
!
!   Write the physical effects specification to the  .res  file.
!   iccutblk() is now an array of length nblock.
!   This is the sixth record on the file
!
      write (imcdf) speed_of_light_au, lfordr, (iccutblk(i), i = 1, symmetry_block_count),            &
                    ltrans, wfact, lvp, lnms, lsms
!
!   Write the grid data to the  .res  file; this is the seventh
!   record on the file
!
      np10 = radial_point_count + 10
      write (imcdf) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, (radial_grid(i), i = 1, np10), &
         & &
                 (radial_grid_derivative(i), i = 1, np10), (radial_derivative_over_radius(i), i = 1, np10) ! (imcdf)
!
!   Write out the interpolated radial wavefunctions; there are
!   2*NW such records; thus the total number of records written
!   at this point is 7+2*NW
!
      do j = 1, orbital_count
         write (imcdf) orbital_energy_parameters(j), orbital_origin_exponents(j), orbital_origin_parameters(j), &
            & orbital_last_grid_index(j)
         write (imcdf) (radial_large_components(i,j), i = 1, orbital_last_grid_index(j)), (radial_small_components(i,j) &
            & , i = 1, orbital_last_grid_index(j))
      enddo
!
      return
      end
