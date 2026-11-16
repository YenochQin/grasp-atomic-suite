!***********************************************************************
!                                                                      *
      subroutine hmout(myid, nprocs, ncf)
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
      use hmat_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: myid
      integer  :: nprocs
      integer, intent(in) :: ncf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ibeg, ico, idiag, list, iro
!-----------------------------------------------
!
!
!  This is original, correct, full matrix
!      ibeg = 1
!      do ico = myid + 1, ncf, nprocs
!         idiag = iendc(ico)
!         do list = ibeg, idiag
!            iro = irow(list)
!            write (99,*) 'H(',iro,ico,')= ', emt(list)
!         enddo
!         ibeg = idiag + 1
!      enddo
!
      ibeg = 1
      do ico = 1, ncf, 4
         ibeg = stored_segment_end_offsets(ico - 1) + 1
         idiag = stored_segment_end_offsets(ico)
         do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            write (96, *) 'H(', iro, ico, ')= ', hamiltonian_values(list)
         end do
         !ibeg = idiag + 1
      end do

      ibeg = stored_segment_end_offsets(1) + 1
      do ico = 2, ncf, 4
         ibeg = stored_segment_end_offsets(ico - 1) + 1
         idiag = stored_segment_end_offsets(ico)
         do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            write (97, *) 'H(', iro, ico, ')= ', hamiltonian_values(list)
         end do
         !ibeg = idiag + 1
      end do

      ibeg = stored_segment_end_offsets(2) + 1
      do ico = 3, ncf, 4
         ibeg = stored_segment_end_offsets(ico - 1) + 1
         idiag = stored_segment_end_offsets(ico)
         do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            write (98, *) 'H(', iro, ico, ')= ', hamiltonian_values(list)
         end do
         !ibeg = idiag + 1
      end do

      ibeg = stored_segment_end_offsets(3) + 1
      do ico = 4, ncf, 4
         ibeg = stored_segment_end_offsets(ico - 1) + 1
         idiag = stored_segment_end_offsets(ico)
         do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            write (99, *) 'H(', iro, ico, ')= ', hamiltonian_values(list)
         end do
         !ibeg = idiag + 1
      end do
      return

      end subroutine hmout
