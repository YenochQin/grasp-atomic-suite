!*******************************************************************
!                                                                  *
      subroutine diaga3(ja1,ja2,ka,ire,iat,rec)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 03  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
!                                                                  *
!   Written by G. Gaigalas,                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use cons_c,          only: zero, one
      use m_c,             only: jlist, jjq1, jjc1, jjc2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ixjtik_i
      use sixj_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ja1, ja2, ka, ire
!      INTEGER, INTENT(OUT) :: IAT
      integer, intent(inout)      :: iat
      real(real64), intent(out) :: rec
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i,ll1,ji,kk1,kk2,iti,iti1,itis,iti1s,ifaz
      real(real64) :: aa, a3
!-----------------------------------------------
      rec = zero
      i=ja1+1
      if(ja1 == 1)i=i+1
      if(ire == 0) then
        if(i >= ja2) then
          rec=one
          iat=1
          return
        end if
        do
          ll1=jlist(i)
          ji=jjq1(3,ll1)-1
          kk2=i-2
          iti=jjc1(kk2)-1
          itis=jjc2(kk2)-1
          kk1=i-1
          iti1=jjc1(kk1)-1
          iti1s=jjc2(kk1)-1
          if(ixjtik(ka,itis,iti,ji,iti1,iti1s) == 0)return
          i=i+1
          if(i == ja2) exit
        end do
        rec=one
        iat=1
        return
      end if
      aa=one
      if(i >= ja2) then
        rec=aa
        iat=1
      else
    1   ll1=jlist(i)
        ji=jjq1(3,ll1)-1
        kk2=i-2
        iti=jjc1(kk2)-1
        itis=jjc2(kk2)-1
        kk1=i-1
        iti1=jjc1(kk1)-1
        iti1s=jjc2(kk1)-1
        if(ire /= 0) then
          call sixj(ka,itis,iti,ji,iti1,iti1s,0,a3)
          a3=a3*sqrt(dble((iti+1)*(iti1s+1)))
          ifaz=ka+ji+iti+iti1s
          if((ifaz/4)*4 /= ifaz)a3=-a3
          aa=aa*a3
        else
          if(ixjtik(ka,itis,iti,ji,iti1,iti1s) == 0)return
        end if
        i=i+1
        if(i == ja2) then
          rec=aa
          iat=1
        else
          go to 1
        end if
      end if
      return
      end subroutine diaga3
