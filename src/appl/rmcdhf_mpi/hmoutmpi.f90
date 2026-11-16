!***********************************************************************
!                                                                      *
      subroutine hmoutmpi(myid, nprocs, ncf)
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
!GG      USE mpi_C
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
      ibeg = 1
      do ico = myid + 1, ncf, nprocs
         idiag = stored_segment_end_offsets(ico)
         do list = ibeg, idiag
            iro = matrix_partner_indices(list)
            write (99,*) 'H(',iro,ico,')= ', hamiltonian_values(list)
         enddo
         ibeg = idiag + 1
      enddo
!
      return
      end subroutine hmoutmpi
