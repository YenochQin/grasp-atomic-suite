!***********************************************************************
!                                                                      *
      subroutine tcsl(n)
!                                                                      *
!   This subroutine transfers data to the initial and final state      *
!   common blocks                                                      *
!                                                                      *
!   Written by Per Jonsson                Last revision: June   1996   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count, orbital_kappa, orbital_l, orbital_principal_n, &
          orbital_symmetry_labels, orbital_two_j
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!
!  Initial state commons
!
!
!  Final state commons
!
!
      if (n == 1) then
         ncfii = csf_count
         nwii = orbital_count
         npii(:orbital_count) = orbital_principal_n(:orbital_count)
         nakii(:orbital_count) = orbital_kappa(:orbital_count)
         nklii(:orbital_count) = orbital_l(:orbital_count)
         nkjii(:orbital_count) = orbital_two_j(:orbital_count)
         nhii(:orbital_count) = orbital_symmetry_labels(:orbital_count)
      else
         ncfff = csf_count
         nwff = orbital_count
         npff(:orbital_count) = orbital_principal_n(:orbital_count)
         nakff(:orbital_count) = orbital_kappa(:orbital_count)
         nklff(:orbital_count) = orbital_l(:orbital_count)
         nkjff(:orbital_count) = orbital_two_j(:orbital_count)
         nhff(:orbital_count) = orbital_symmetry_labels(:orbital_count)
      end if

      return
      end subroutine tcsl
