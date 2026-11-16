!***********************************************************************
      subroutine spicmvmpi (n,m,b,c)
!     IMPLICIT REAL*8          (A-H,O-Z)

!   This routine now works for both rscfmpivu and rcimpivu. By removing
!   the include statement and the call to gdsummpi and setting myid=0,
!   nprocs=1, it also works for the corresponding serial program.
!   Matrix is stored in the mode of upper-triangle-by columns, or
!   you can say lower-triangle-by-rows. 98-08-06
!
!   Matrix-matrix product: C = AB.  A  sparse  representation of the   *
!   lower triangle of the  (NxN)  matrix  A  is assumed available in   *
!   COMMON/HMAT/.                                                      *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos'  routine  SPSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT                                 *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 19 Dec 1992   *
!   MPI version by Xinghong He            Last revision: 29 Jul 1998   *
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use mpi_c
      use hmat_c,          only: hamiltonian_values, stored_segment_end_offsets, matrix_partner_indices, stored_matrix_element_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dinit_i
!      USE dmerge_I
      use dmerge_i, only: dmerge
      use mpiu_i, only: gdsummpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(real64), intent(in) :: b(n,m)
      real(real64) :: c(n,m)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icol, ibeg, iend, nelc, iv
      real(real64) :: diag, dl
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!
!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- no overdimensioning
!
      call dinit (n*m, 0.0_real64, c, 1)

      ibeg = 1
      do icol = mpi_rank + 1, n, mpi_process_count
            !IBEG = IENDC(ICOL-1)+1
            !IEND = IENDC(ICOL)
         iend = stored_segment_end_offsets(icol)
         nelc = iend - ibeg + 1
         do iv = 1, m
            diag =  c(icol,iv) + hamiltonian_values(iend)*b(icol,iv)
            call dmerge (nelc-1,b(:n,iv),c(:n,iv),               &
                         matrix_partner_indices(ibeg:iend-1),hamiltonian_values(ibeg:iend-1),b(icol,iv),dl)
            c(icol,iv) = diag + dl
         enddo
         ibeg = iend + 1
      enddo

      call gdsummpi (c, n*m)

      return
      end subroutine spicmvmpi
