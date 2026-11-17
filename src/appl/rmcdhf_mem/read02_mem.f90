!*******************************************************************
!                                                                  *
      subroutine read02_mem(jblock)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use hmat_c,     only: matrix_partner_indices
      use rmcdhf_mem_c, only: imax2_30, imin2_30, irow_30
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jblock
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------

      do i = imin2_30(jblock)+1, imax2_30(jblock)
         matrix_partner_indices(i-imin2_30(jblock)) = irow_30(i)
      end do

      return
      end subroutine read02_mem
