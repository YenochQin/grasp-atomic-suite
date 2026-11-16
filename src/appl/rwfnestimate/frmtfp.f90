

!***********************************************************************
!                                                                      *
      subroutine frmtfp(index, nsubs)
!                                                                      *
!   This  subroutine  is  used  to  produce  estimates  of  the wave   *
!   functions by use of the Thomas-Fermi approximation to the direct   *
!   potential. SUBROUTINE  SOLVH  is used  to obtain the radial wave   *
!   functions.                                                         *
!                                                                      *
!   Call(s) to: [ERWF]: SOLVH.                                         *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:06:21   1/ 2/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use parameter_def, only: max_radial_points, radial_workspace_capacity
      use vast_kind_param, only:  double
      use left_c, only: set
      use orb_c, only: orbital_principal_n, orbital_symmetry_labels
      use wave_c
      use whfrom_c, only: source
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use solvh_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: nsubs
      integer , intent(in) :: index(max_orbitals)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, loc
      logical :: fail
!-----------------------------------------------
!
!
      do j = 1, nsubs
!
         loc = index(j)
!
         if (set(loc)) cycle
!
!   Estimate the leading coefficient of the series expansion for
!   the orbital
!
         orbital_origin_parameters(loc) = 10.0d00
!
!   Calculate radial wavefunctions
!
         call solvh (loc, fail)
!
!   Message if SOLVH did not converge; reset SET(LOC) otherwise
!
         if (fail) then
            write (*, 300) orbital_principal_n(loc), orbital_symmetry_labels(loc)
         else
            set(loc) = .true.
               !SOURCE(LOC) = 'Thomas-Fermi estimate'
            source(loc) = 'T-F'
         endif
!
      end do
!
      return
!
  300 format(/,'TFWAVE: Unable to compute radial'/,/,' wavefunction for ',i2,a2&
         ,' subshell;')
      return
!
      end subroutine frmtfp
