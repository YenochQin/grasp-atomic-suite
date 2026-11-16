!***********************************************************************
!                                                                      *
      subroutine connect
!                                                                      *
!   The position of an orbital in the merged list is connected to      *
!   the positions in the initial and final state lists                 *
!                                                                      *
!   Written by Per Jonsson                                             *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use biorb_c
      use orb_c
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
!-----------------------------------------------
!
!ww      INTEGER PNTRIQ
!     POINTER (PNTRIQ,RIQDUMMY)
!
!
!     Initialize
!
      nnii(:orbital_count) = 0.0_real64
      nnff(:orbital_count) = 0.0_real64
!
!   Loop over the orbitals in the merged list
!
      do i = 1, orbital_count
         do j = 1, nwii
            if (orbital_principal_n(i)/=npii(j) .or. orbital_kappa(i)/=nakii(j)) cycle
            nnii(i) = j
         end do
         do j = 1, nwff
            if (orbital_principal_n(i)/=npff(j) .or. orbital_kappa(i)/=nakff(j)) cycle
            nnff(i) = j
         end do
      end do

      return
      end subroutine connect
