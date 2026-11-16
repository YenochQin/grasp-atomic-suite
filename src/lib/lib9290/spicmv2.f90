!***********************************************************************
      subroutine spicmv2(n, m, b, c)
!
!  Modified from the mpi version spicmvmpi.f by simply removing things
!  required by mpi communication.
!
!  Xinghong He  98-07-29
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:47   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use hmat_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dinit_i
      use dmerge_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(double), dimension(n,m), intent(in) :: b
      real(double), dimension(n,m), intent(inout) :: c
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ibeg, icol, iend, nelc, iv
      real(double) :: diag, dl
!-----------------------------------------------
!
!
!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- no overdimensioning
!
      call dinit (n*m, 0.d0, c, 1)

      ibeg = 1
      do icol = 1, n
            !IBEG = IENDC(ICOL-1)+1
            !IEND = IENDC(ICOL)
         iend = stored_segment_end_offsets(icol)
         nelc = iend - ibeg + 1
         do iv = 1, m
            diag = c(icol,iv) + hamiltonian_values(iend)*b(icol,iv)
            call dmerge (nelc-1, b(:n,iv), c(:n,iv), matrix_partner_indices(ibeg:iend-1), hamiltonian_values(ibeg:iend-1), b(&
               icol,iv), dl)
            c(icol,iv) = diag + dl
         end do
         ibeg = iend + 1
      end do

      return
      end subroutine spicmv2
