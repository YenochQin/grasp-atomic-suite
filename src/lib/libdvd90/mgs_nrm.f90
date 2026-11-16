module mgs_nrm_i
    implicit none (type, external)
    private
    public :: mgs_nrm
contains
      subroutine mgs_nrm(n, kp, new, scra, b)
!=======================================================================
!       Orthogonalizis the new vectors in B (from kp+1...kp+new)
!     to the previous kp B vectors and to themselves
!     using Modified Gram Schmidt. Then normalizes them.
!       The procedure is repeated twice.
!-----------------------------------------------------------------------
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
!     USE mpi_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      !use dgemv_I
      !use ddot_I
      !use dscal_I
      use suite_linear_algebra, only: ddot, dgemv, dscal
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer , intent(in) :: kp
      integer, intent(in) :: new
      real(real64), intent(inout) :: scra(new)
      real(real64), intent(inout) :: b((kp + new)*n)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, newstart, kcur, k, jcur, j, mm
      real(real64) :: dnm
!-----------------------------------------------

!
! MODIFIED GRAM-SCHMIDT  (twice)
!
        !if (myid.eq.0) then
      do i = 1, 2

         newstart = kp*n + 1
!
!  First record contribution from the kp vectors to the new ones.
!
         kcur = 1
         do k = 1, kp
            jcur = newstart
            call dgemv ("T", n, new, 1.0_real64, b(jcur), n, b(kcur), 1, 0.0_real64, scra, &
               1)
            do j = 1, new
!           call daxpy(N,-scra(j),B(kcur),1,B(jcur),1)
               do mm = 0, n - 1
                  b(jcur+mm) = b(jcur+mm) - scra(j)*b(kcur+mm)
               end do
               jcur = jcur + n
            end do
            kcur = kcur + n
         end do
!
!  Then orthogonalize the new ones among themselves.
!
         do k = 1, new
            jcur = kcur + n
!  The current vector should be normalized
!
            dnm = ddot(n,b(kcur),1,b(kcur),1)
            dnm = sqrt(dnm)
            call dscal (n, 1/dnm, b(kcur), 1)
!
            ! No trailing vectors: avoid forming B one element past its end.
            if (k < new) then
            call dgemv ("T", n, new - k, 1.0_real64, b(jcur), n, b(kcur), 1, 0.0_real64, &
               scra, 1)
            end if
            do j = k + 1, new
!              call daxpy(N,-scra(j-k),B(kcur),1,B(jcur),1)
               do mm = 0, n - 1
                  b(jcur+mm) = b(jcur+mm) - scra(j-k)*b(kcur+mm)
               end do
               jcur = jcur + n
            end do
            kcur = kcur + n
         end do
      end do
        !end if

!        call MPI_BCAST(B,(kp+New)*n,MPI_DOUBLE_PRECISION,0,
!     :                 MPI_COMM_WORLD,ierr)

      return
      end subroutine mgs_nrm
end module mgs_nrm_i
