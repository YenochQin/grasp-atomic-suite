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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use hmat_c,          only: hamiltonian_values
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dinit_i
!      USE dmerge_dnicmv_I
      use dmerge_dnicmv_i, only: dmerge_dnicmv
      use mpiu_i, only: gdsummpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(real64), intent(in) :: b(n, m)
      real(real64) :: c(n, m)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: ibeg, iend, icol, nelc, iv
      real(real64) :: diag, dl
!-----------------------------------------------

!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- there is no overdimensioning

      call dinit (n*m, 0.0_real64, c, 1)

      ibeg = 1
      iend = 0
      do icol = mpi_rank + 1, n, mpi_process_count
         iend = iend + icol
         nelc = iend - ibeg + 1
         do iv = 1, m
            diag = c(icol,iv) + hamiltonian_values(iend)*b(icol,iv)
            call dmerge_dnicmv (nelc - 1, b(1:n,iv), c(1:n,iv), &
                                hamiltonian_values(ibeg:iend), b(icol,iv), dl)
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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      real(real64), intent(in) :: dconst
      real(real64), intent(out) :: dl
      real(real64), intent(in) :: db(*)
      real(real64), intent(inout) :: dc(*)
      real(real64), intent(in) :: da(n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
      real(real64) :: dsum
!-----------------------------------------------

      dsum = 0._real64
      dsum = dot_product(da,db(:n))
      dc(:n) = dc(:n) + dconst*da
      dl = dsum

      return
      end subroutine dmerge_dnicmv
