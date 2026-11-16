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
      use vast_kind_param, only:  double
      use mpi_c
      use hmat_c,          only: emt, iendc, irow, nelmnt
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dinit_i
!      USE dmerge_I
      use dmerge_i, only: dmerge
      use mpiu_i, only: gdsummpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: n
      integer, intent(in) :: m
      real(double) :: b(n,m)
      real(double) :: c(n,m)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: icol, ibeg, iend, nelc, iv
      real(double) :: diag, dl
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!
!   Initialise the result matrix; note that this is specific to the
!   data structure of DVDSON --- no overdimensioning
!
      call dinit (n*m, 0.0d00, c, 1)

      ibeg = 1
      do icol = myid + 1, n, nprocs
            !IBEG = IENDC(ICOL-1)+1
            !IEND = IENDC(ICOL)
         iend = iendc(icol)
         nelc = iend - ibeg + 1
         do iv = 1, m
            diag =  c(icol,iv) + emt(iend)*b(icol,iv)
            call dmerge (nelc-1,b(:n,iv),c(:n,iv),               &
                         irow(ibeg:iend-1),emt(ibeg:iend-1),b(icol,iv),dl)
            c(icol,iv) = diag + dl
         enddo
         ibeg = iend + 1
      enddo

      call gdsummpi (c, n*m)

      return
      end subroutine spicmvmpi
