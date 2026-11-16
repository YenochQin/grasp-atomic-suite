!*******************************************************************
!                                                                  *
      subroutine read03_mem(jblock)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use hmat_c,     only: irow, iendc
      use rmcdhf_mem_c
      use mpi_s
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: jblock
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------

      do i = imin1_30(jblock)+1+myid, imax1_30(jblock), nprocs
         iendc(i-imin1_30(jblock)) = iendc_30(i)
      end do
      do i = imin2_30(jblock)+1, imax2_30(jblock)
         irow(i-imin2_30(jblock)) = irow_30(i)
      end do

      return
      end subroutine read03_mem
