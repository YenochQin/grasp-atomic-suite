module gracah1_i
    implicit none (type, external)
    private
    public :: gracah1
contains
!*******************************************************************
!                                                                  *
      subroutine gracah1(i,j,k,l,m,n,rac)
!                                                                  *
!   SUBROUTINE TO CALCULATE RACAH COEFFICIENTS.                    *
!   THE ARGUMENTS I,J,K,L,M,N SHOULD BE TWICE THEIR ACTUAL VALUE.  *
!   WRITTEN BY N. S. SCOTT                                         *
!   Modified by C. Froese Fischer, March 11, 1988 to use table     *
!                                                         look-up  *
!   and                                                            *
!   G. Gaigalas, September 14, 1997                                *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, one, two
      use facts_c
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: i
      integer , intent(in) :: j
      integer , intent(in) :: k
      integer , intent(in) :: l
      integer , intent(in) :: m
      integer , intent(in) :: n
      real(real64) , intent(out) :: rac
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j1,j2,j3,j4,j5,j6,j7,numin,numax,icount,kk,ki
!-----------------------------------------------
      j1 = i+j+m
      j2 = k+l+m
      j3 = i+k+n
      j4 = j+l+n
      if (mod(j1,2) == 0  .and.  mod(j2,2) == 0   .and. &
          mod(j3,2) == 0  .and.  mod(j4,2) == 0 )  then
          j1 = j1/2
          j2 = j2/2
          j3 = j3/2
          j4 = j4/2
          if (max(i,j,m) <= j1 .and.  max(k,l,m) <= j2  .and. &
              max(i,k,n) <= j3 .and.  max(j,l,n) <= j4  )  then
              j5 = (i+j+k+l)/2
              j6 = (i+l+m+n)/2
              j7 = (j+k+m+n)/2
              numin = max(j1, j2, j3, j4) + 1
              numax = min(j5, j6, j7)     + 1
              rac = one
              icount = 0
              do kk = numin+1,numax
                ki = numax - icount
                rac = one - (rac*(ki*(j5-ki+2)*(j6-ki+2)*(j7-ki+2)))/ &
                         ((ki-1-j1)*(ki-1-j2)*(ki-1-j3)*(ki-1-j4))
                icount = icount+1
              end do
              rac = rac*exp(                                          &
                    (gam(numin+1) - gam(numin-j1) - gam(numin-j2) -   &
                     gam(numin-j3) - gam(numin-j4) - gam(j5+2-numin)- &
                     gam(j6+2-numin)-gam(j7+2-numin)) +               &
                    (gam(j1+1-i)+gam(j1+1-j)+gam(j1+1-m)-gam(j1+2) +  &
                     gam(j2+1-k)+gam(j2+1-l)+gam(j2+1-m)-gam(j2+2) +  &
                     gam(j3+1-i)+gam(j3+1-k)+gam(j3+1-n)-gam(j3+2) +  &
                     gam(j4+1-j)+gam(j4+1-l)+gam(j4+1-n)-gam(j4+2))/two)
              if (mod(j5+numin,2) == 0) rac = -rac
          else
             rac = zero
          end if
      else
         rac = zero
      end if
      return
      end
end module gracah1_i
