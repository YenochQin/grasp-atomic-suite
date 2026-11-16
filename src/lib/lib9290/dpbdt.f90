!***********************************************************************
!                                                                      *
      subroutine dpbdt(j)
!-----------------------------------------------
!                                                                      *
!   This subroutine computes H times the derivative, with respect to   *
!   the internal grid, of the large and small components of the wave   *
!   function with index  J .  These  are tabulated, respectively, in   *
!   arrays  TA  and  TB  in  COMMON  block  /TATB/ .                   *
!                                                                      *
!   A  thirteen-point  Lagrange formaula is used for the calculation   *
!   of derivatives.                                                    *
!                                                                      *
!   Written by Farid F Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:47:21   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use grid_c
      use lic13_c
      use tatb_c, only: ta, tb
      use wave_c
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer  :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, irow, loc
      real(double) :: a1, a2, a3, a4, a5, a6, hdpbdt, hdqbdt, aik
!-----------------------------------------------
!

!     Replace the equivalence statement
      a1=a(7,1)
      a2=a(7,2)
      a3=a(7,3)
      a4=a(7,4)
      a5=a(7,5)
      a6=a(7,6)
!
!   Compute derivative in three separate regions
!
!   First, points 1 to 6
!
      do i = 1,6
         hdpbdt = 0.d0
         hdqbdt = 0.d0
!cjb replace loop by SUM
!        DO 1 K = 1,13
!            AIK = A(I,K)
!            HDPBDT = HDPBDT+AIK*PF(K,J)
!            HDQBDT = HDQBDT+AIK*QF(K,J)
!    1    CONTINUE
         hdpbdt = sum(a(i   ,:)*pf(:13,j))
         hdqbdt = sum(a(i   ,:)*qf(:13,j))
!cjb
         ta(i) = hdpbdt
         tb(i) = hdqbdt
      end do

!   Next, points 7 to N-6
!
!   Special treatment for this region because of the symmetry of
!   the differentiation formula
!
      do i = 7, n - 6
         ta(i) = a1*(pf(i - 6,j) - pf(i + 6,j)) + a2*(pf(i - 5,j) - pf(i + 5,j)&
            ) + a3*(pf(i - 4,j) - pf(i + 4,j)) + a4*(pf(i - 3,j) - pf(i + 3,j))&
             + a5*(pf(i - 2,j) - pf(i + 2,j)) + a6*(pf(i - 1,j) - pf(i + 1,j))
         tb(i) = a1*(qf(i - 6,j) - qf(i + 6,j)) + a2*(qf(i - 5,j) - qf(i + 5,j)&
            ) + a3*(qf(i - 4,j) - qf(i + 4,j)) + a4*(qf(i - 3,j) - qf(i + 3,j))&
             + a5*(qf(i - 2,j) - qf(i + 2,j)) + a6*(qf(i - 1,j) - qf(i + 1,j))
      end do
!
!   Last, points N-5 to N
!
      do i = n - 5, n
         irow = i - n + 13
          hdpbdt = sum(a(irow,:)*pf(n-12:n,j))
          hdqbdt = sum(a(irow,:)*qf(n-12:n,j))
!- git 1/12/07 code below from the new libraries
!         HDPBDT = 0.D0
!        HDQBDT = 0.D0
!        DO K = 1, 13
!           AIK = A(IROW,K)
!           LOC = N - 13 + K
!           HDPBDT = HDPBDT + AIK*PF(LOC,J)
!           HDQBDT = HDQBDT + AIK*QF(LOC,J)
!        END DO
!- git 1/12/07
         ta(i) = hdpbdt
         tb(i) = hdqbdt
      end do
      return
      end subroutine dpbdt
