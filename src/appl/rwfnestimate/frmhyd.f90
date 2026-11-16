

!***********************************************************************
!                                                                      *
      subroutine frmhyd(index, nsubs, modify, nradial)
!                                                                      *
!   This  subroutine  is  used  to  produce  estimates  of  the wave   *
!   functions from the hydrogenic approximation. The screening cons-   *
!   tant is taken to be 0.4 (Changed, see comments below,XHH)          *
!                                                                      *
!   Call(s) to: [ERWF]: DCBSRW.                                        *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: nnn1
      use vast_kind_param, only:  double
      use def_c
      use left_c
      use orb_c, only: e, np, nak, nh
      use wave_c, only: pf, qf, pz, mf
      use whfrom_c, only: source
      use hydpar_c, only: sigma
      use iounit_c
      use default_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dcbsrw_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nsubs
      logical , intent(in) :: modify
      integer , intent(in) :: index(nnnw)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, loc, nradial
      real(double) :: zeff, delz
      logical :: yes, getyn
!-----------------------------------------------
!
!
      delz = 0.d0
      if (nradial == 4)  then
         write (istde, *) 'Enter increase in Z for correlation orbitals'
         read  (5, *) delz
      end if
      write (istde, *) 'Orbital Z_eff for hydrogenic orbitals'
      do j = 1, nsubs
         loc = index(j)
         if (set(loc)) cycle
         if (modify) then
            write (istde, *) 'Input new value >'
            read (5, *) sigma(loc)
         endif
         zeff = z - sigma(loc) + delz
         write (istde, '(I2,A2,F10.2)') np(loc), nh(loc), zeff
!           ...Calculate radial wavefunctions
         call dcbsrw (np(loc), nak(loc), zeff, e(loc), pz(loc), pf(:,loc), &
            qf(:,loc), mf(loc))
         set(loc) = .true.
         source(loc) = 'Hyd'
      end do
      write (istde, *)
      return
      end subroutine frmhyd
