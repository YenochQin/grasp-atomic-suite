

!***********************************************************************
!                                                                      *
      subroutine frmrwf(index, nsubs, filnam)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, OPENFL.                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use debug_c
      use grid_c
      use left_c
      use orb_c, only: orbital_energy_parameters, orbital_principal_n, orbital_kappa, orbital_symmetry_labels
      use wave_c, only: orbital_origin_parameters
      use whfrom_c, only: source
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use intrpq_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nsubs
      character  :: filnam*(*)
      integer , intent(in) :: index(max_orbitals)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr, ios, npy, naky, my, j, loc, i, lenth
      real(real64) :: ey, dnorm
      real(real64), pointer :: pa(:), qa(:), ra(:)
      logical :: found
      character :: g92rwf*6
!-----------------------------------------------
!
!

      !FORM = 'UNFORMATTED'
      !STATUS = 'OLD'

      call openfl (23, filnam, 'UNFORMATTED', 'OLD', ierr)
      if (ierr == 1) then
         write (stderr_unit, *) 'Error openning file "', filnam(1:len_trim(filnam)), &
            '"'
         close(23)
         stop
      endif
!
!   Check the file; if not as expected, try again
!
      read (23, iostat=ios) g92rwf
      if (ios/=0 .or. g92rwf/='G92RWF') then
         write (stderr_unit, *) 'This is not a Radial WaveFunction File;'
         close(23)
         stop
      endif
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      if (ldbpr(3)) write (99, 300)
    2 continue
      found = .false.
      read (23, iostat=ios) npy, naky, ey, my
      if (ios == 0) then
         do j = 1, nsubs
            loc = index(j)
            if (.not.(.not.set(loc) .and. orbital_principal_n(loc)==npy .and. orbital_kappa(loc)==naky)) &
               cycle
            found = .true.
            orbital_energy_parameters(loc) = ey
            call alloc (pa, my, 'PA', 'FRMFRW')
            call alloc (qa, my, 'QA', 'FRMFRW')
            call alloc (ra, my, 'RA', 'FRMFRW')
            read (23) orbital_origin_parameters(loc), (pa(i),i=1,my), (qa(i),i=1,my)
            read (23) (ra(i),i=1,my)
            call intrpq (pa, qa, my, ra, loc, dnorm)
            if (ldbpr(3)) write (99, 301) orbital_principal_n(loc), orbital_symmetry_labels(loc), &
               & orbital_energy_parameters(loc), dnorm
            call dalloc (pa, 'PA', 'FRMFRW')
            call dalloc (qa, 'QA', 'FRMFRW')
            call dalloc (ra, 'RA', 'FRMFRW')
            lenth = len_trim(filnam)
            set(loc) = .true.
            source(loc) = filnam(1:3)
            go to 2
         end do
         if (.not.found) then
            read (23)
            read (23)
            go to 2
         endif
      endif
      if (ldbpr(3)) write (99, *) ' orbitals renormalised;'
!
      close(23)
!
      return
!
  300 format(/,'From SUBROUTINE FRMRWF:'/,' Orbital',8x,'Eigenvalue',19x,'Norm'&
         )
  301 format(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)
      return
!
      end subroutine frmrwf
