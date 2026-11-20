

!***********************************************************************
!                                                                      *
      subroutine frmtfp(matched_orbital_indices, matched_subshell_count)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use grasp_parameters, only: max_radial_points, radial_workspace_capacity
      use, intrinsic :: iso_fortran_env, only: real64
      use initial_orbital_selection, only: set
      use orbital_csf_basis, only: orbital_principal_n, orbital_symmetry_labels
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use initial_orbital_source, only: source
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use solvh_i, only: solvh
      implicit none (type, external)
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: matched_subshell_count
      integer , intent(in) :: matched_orbital_indices(max_orbitals)
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
      do j = 1, matched_subshell_count
!
         loc = matched_orbital_indices(j)
!
         if (set(loc)) cycle
!
!   Estimate the leading coefficient of the series expansion for
!   the orbital
!
         orbital_origin_parameters(loc) = 10.0_real64
!
!   Calculate radial wavefunctions
!
         call solvh (loc, fail)
!
!   Message if SOLVH did not converge; reset SET(LOC) otherwise
!
         if (fail) then
            write (*, '(/,"TFWAVE: Unable to compute radial"/,/," wavefunction for ",i2,a2," subshell;")') &
               orbital_principal_n(loc), orbital_symmetry_labels(loc)
         else
            set(loc) = .true.
               !SOURCE(LOC) = 'Thomas-Fermi estimate'
            source(loc) = "T-F"
         end if
!
      end do
!
      return
!
      return
!
      end subroutine frmtfp
