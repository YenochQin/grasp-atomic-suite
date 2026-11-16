!***********************************************************************
!                                                                      *
      subroutine dnicmv(n, m, b, c)
!                                                                      *
!   Matrix-matrix product: C = AB.  The lower triangle of the  (NxN)   *
!   matrix is assumed available in packed form in the array EMT. The   *
!   matrices B and C are (NxM).                                        *
!                                                                      *
!   This is an adaptation of  Andreas Stathopulos'  routine  TRSBMV,   *
!   and is specific to GRASP2 derivatives.                             *
!                                                                      *
!   Call(s) to: [AUXBLAS]: DINIT/SINIT;                                *
!               [BLAS]: DAXPY/SAXPY, DDOT/SDOT.                        *
!                                                                      *
!   F A Parpia and A Stathopoulos         Last revision: 09 Oct 1992   *
!   Block version by Xinghong He          Last revision: 18 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use hmat_c,          only: emt
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dinit_i
!      USE dmerge_dnicmv_I
      use dmerge_dnicmv_i, only: dmerge_dnicmv
      use mpiu_i, only: gdsummpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(double), dimension(n,m) :: b
      real(double), dimension(n,m) :: c
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ibeg, iend, icol, nelc, iv
      real(double) :: diag, dl
!-----------------------------------------------

!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- there is no overdimensioning

      call dinit (n*m, 0.0d00, c, 1)

      ibeg = 1
      iend = 0
      do icol = myid + 1, n, nprocs
         iend = iend + icol
         nelc = iend - ibeg + 1
         do iv = 1, m
            diag = c(icol,iv) + emt(iend)*b(icol,iv)
            call dmerge_dnicmv (nelc - 1, b(1:n,iv), c(1:n,iv), &
                                emt(ibeg:iend), b(icol,iv), dl)
            c(icol,iv) = diag + dl
         end do
         ibeg = iend + 1
      end do

      call gdsummpi (c,n*m)

      return
      end subroutine dnicmv



!***********************************************************************
!                                                                      *
      subroutine dmerge_dnicmv(n, db, dc, da, dconst, dl)
!
!  Used by dnimcv
!  Developed from dmerge. The only diff is: idy not needed here
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      real(double), intent(in) :: dconst
      real(double), intent(out) :: dl
      real(double), dimension(*), intent(in) :: db
      real(double), dimension(*), intent(inout) :: dc
      real(double), dimension(n), intent(in) :: da
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(double) :: dsum
!-----------------------------------------------

      dsum = 0.d0
      dsum = dot_product(da,db(:n))
      dc(:n) = dc(:n) + dconst*da
      dl = dsum

      return
      end subroutine dmerge_dnicmv
