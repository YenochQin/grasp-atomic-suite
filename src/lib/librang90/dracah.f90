module dracah_i
    implicit none (type, external)
    private
    public :: dracah
contains
!***********************************************************************
!                                                                      *
      subroutine dracah(i, j, k, l, m, n, rac)
!                                                                      *
!   SUBROUTINE  to calculate Racah coefficients. The arguments I, J,   *
!   K, L, M, N should be twice their actual value. Works for integer   *
!   and  half-integer  values of  angular momenta. The routine makes   *
!   use of the GAM  array, thus  SUBROUTINE FACTT must be called be-   *
!   fore this routine is used.                                         *
!                                                                      *
!   Written by N S Scott                    Last update: 16 Oct 1992   *
!   The last modification made by G. Gaigalas           October 2017   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use facts_c, only: double, gam, mfact
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
      integer :: j1, j2, j3, j4, j5, j6, j7, numin, numax, icount, kk, ki
!-----------------------------------------------
      j1 = i + j + m
      j2 = k + l + m
      j3 = i + k + n
      j4 = j + l + n
      if (2*max(max(i,j),m) - j1>0 .or. mod(j1,2)/=0) goto 2
      if (2*max(max(k,l),m) - j2>0 .or. mod(j2,2)/=0) goto 2
      if (2*max(max(i,k),n) - j3>0 .or. mod(j3,2)/=0) goto 2
      if (2*max(max(j,l),n) - j4>0 .or. mod(j4,2)/=0) goto 2
      goto 1
    2 continue
      rac = 0.0_real64
      return
!
    1 continue
      j1 = j1/2
      j2 = j2/2
      j3 = j3/2
      j4 = j4/2
      j5 = (i + j + k + l)/2
      j6 = (i + l + m + n)/2
      j7 = (j + k + m + n)/2
      numin = max(max(max(j1,j2),j3),j4) + 1
      numax = min(min(j5,j6),j7) + 1
      rac = 1.0_real64
      icount = 0
!
      if (numin /= numax) then
         numin = numin + 1
!
         do kk = numin, numax
            ki = numax - icount
            rac = 1.0_real64 - rac*dble(ki*(j5 - ki + 2)*(j6 - ki + 2)*(j7 - ki + 2&
               ))/dble((ki - 1 - j1)*(ki - 1 - j2)*(ki - 1 - j3)*(ki - 1 - j4))
            icount = icount + 1
         end do
!
         numin = numin - 1
      end if
      rac = rac*(-1.0_real64)**(j5 + numin + 1)*exp((gam(numin+1)-gam(numin-j1)-gam&
         (numin-j2)-gam(numin-j3)-gam(numin-j4)-gam(j5+2-numin)-gam(j6+2-numin)&
         -gam(j7+2-numin))+(gam(j1+1-i)+gam(j1+1-j)+gam(j1+1-m)-gam(j1+2)+gam(&
         j2+1-k)+gam(j2+1-l)+gam(j2+1-m)-gam(j2+2)+gam(j3+1-i)+gam(j3+1-k)+gam(&
         j3+1-n)-gam(j3+2)+gam(j4+1-j)+gam(j4+1-l)+gam(j4+1-n)-gam(j4+2))*&
         0.5_real64)
!
      return
      end subroutine dracah
end module dracah_i
