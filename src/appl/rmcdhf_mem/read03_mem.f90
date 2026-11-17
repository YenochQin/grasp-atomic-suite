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
      use hmat_c,     only: matrix_partner_indices, stored_segment_end_offsets
      use rmcdhf_mem_c, only: iendc_30, imax1_30, imax2_30, imin1_30, imin2_30, irow_30
      use mpi_s, only: mpi_process_count, mpi_rank
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

      do i = imin1_30(jblock)+1+mpi_rank, imax1_30(jblock), mpi_process_count
         stored_segment_end_offsets(i-imin1_30(jblock)) = iendc_30(i)
      end do
      do i = imin2_30(jblock)+1, imax2_30(jblock)
         matrix_partner_indices(i-imin2_30(jblock)) = irow_30(i)
      end do

      return
      end subroutine read03_mem
