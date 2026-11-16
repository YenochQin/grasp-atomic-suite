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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:08:49   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use biorb_c
      use orb_c
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
      endif

      return
      end subroutine tcsl
