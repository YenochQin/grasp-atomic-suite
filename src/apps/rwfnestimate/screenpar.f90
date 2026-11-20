

!***********************************************************************
      subroutine screenpar(core_subshell_count)
!
! Purpose:
!   Compute hydrogenic screen parameters
!
! Input:
!   ncore
!   in the common - nw, nkj()
!
! Output:
!   in the common - sigma()
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis, only: orbital_count, orbital_two_j
      use hydrogenic_screening, only: sigma
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: core_subshell_count
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nelectron, i
!-----------------------------------------------

       !...Core orbitals
      nelectron = 0
      do i = 1, core_subshell_count
         sigma(i) = nelectron + (orbital_two_j(i)+1)/2
         nelectron = nelectron + orbital_two_j(i) + 1
      end do

       !...Peel orbitals
      sigma(core_subshell_count+1:orbital_count) = nelectron

      return
      end subroutine screenpar
