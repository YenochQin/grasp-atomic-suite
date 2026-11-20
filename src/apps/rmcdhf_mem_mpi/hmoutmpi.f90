!***********************************************************************
!                                                                      *
      subroutine hmoutmpi(mpi_rank, mpi_process_count, block_csf_count)
!                                                                      *
!   Routine for printing the Hamiltonian matrix.                       *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 21 Dec 1992   *
!   Block Version by Xinghong He          Last revision: 30 Jan 1999   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
!GG      USE mpi_C
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: mpi_rank
      integer  :: mpi_process_count
      integer, intent(in) :: block_csf_count
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ibeg, ico, idiag, list, iro
!-----------------------------------------------
!
      ibeg = 1
      do ico = mpi_rank + 1, block_csf_count, mpi_process_count
         idiag = stored_segment_end_offsets(ico)
         do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            write (99,*) "H(",iro,ico,")= ", hamiltonian_values(list)
         end do
         ibeg = idiag + 1
      end do
!
      return
      end subroutine hmoutmpi
