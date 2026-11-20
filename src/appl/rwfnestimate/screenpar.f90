

!***********************************************************************
      subroutine screenpar(ncore)
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
      integer , intent(in) :: ncore
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nelectron, i
!-----------------------------------------------

       !...Core orbitals
      nelectron = 0
      do i = 1, ncore
         sigma(i) = nelectron + (orbital_two_j(i)+1)/2
         nelectron = nelectron + orbital_two_j(i) + 1
      end do

       !...Peel orbitals
      sigma(ncore+1:orbital_count) = nelectron

      return
      end subroutine screenpar
