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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grid_c
      use lic13_c
      use tatb_c, only: ta, tb
      use wave_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, k, irow, loc
      real(real64) :: a1, a2, a3, a4, a5, a6, hdpbdt, hdqbdt, aik
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
         hdpbdt = 0._real64
         hdqbdt = 0._real64
!cjb replace loop by SUM
!        DO 1 K = 1,13
!            AIK = A(I,K)
!            HDPBDT = HDPBDT+AIK*PF(K,J)
!            HDQBDT = HDQBDT+AIK*QF(K,J)
!    1    CONTINUE
         hdpbdt = sum(a(i   ,:)*radial_large_components(:13,j))
         hdqbdt = sum(a(i   ,:)*radial_small_components(:13,j))
!cjb
         ta(i) = hdpbdt
         tb(i) = hdqbdt
      end do

!   Next, points 7 to N-6
!
!   Special treatment for this region because of the symmetry of
!   the differentiation formula
!
      do i = 7, radial_point_count - 6
         ta(i) = a1*(radial_large_components(i - 6,j) - radial_large_components(i + 6,j)) + a2* &
            & (radial_large_components(i - 5,j) - radial_large_components(i + 5,j)&
            ) + a3*(radial_large_components(i - 4,j) - radial_large_components(i + 4,j)) + a4* &
               & (radial_large_components(i - 3,j) - radial_large_components(i + 3,j))&
             + a5*(radial_large_components(i - 2,j) - radial_large_components(i + 2,j)) + a6* &
                & (radial_large_components(i - 1,j) - radial_large_components(i + 1,j))
         tb(i) = a1*(radial_small_components(i - 6,j) - radial_small_components(i + 6,j)) + a2* &
            & (radial_small_components(i - 5,j) - radial_small_components(i + 5,j)&
            ) + a3*(radial_small_components(i - 4,j) - radial_small_components(i + 4,j)) + a4* &
               & (radial_small_components(i - 3,j) - radial_small_components(i + 3,j))&
             + a5*(radial_small_components(i - 2,j) - radial_small_components(i + 2,j)) + a6* &
                & (radial_small_components(i - 1,j) - radial_small_components(i + 1,j))
      end do
!
!   Last, points N-5 to N
!
      do i = radial_point_count - 5, radial_point_count
         irow = i - radial_point_count + 13
          hdpbdt = sum(a(irow,:)*radial_large_components(radial_point_count-12:radial_point_count,j))
          hdqbdt = sum(a(irow,:)*radial_small_components(radial_point_count-12:radial_point_count,j))
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
