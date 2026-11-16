module diaga4_i
    implicit none (type, external)
    private
    public :: diaga4
contains
!*******************************************************************
!                                                                  *
      subroutine diaga4(ja1,ja2,k1,k2,ka,ire,iat,recc)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 04  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  NINE                                     *
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
      use cons_c,          only: zero
      use m_c,             only: jlist, jjq1, jjq2, jjc1, jjc2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i
      use nine_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ja1,ja2,k1,k2,ka,ire
!      INTEGER, INTENT(OUT)      :: IAT
      integer, intent(inout)      :: iat
      real(real64), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ij1,ij2,ia1,ia2,ib1,ib2,n1,n2,j2,j2s,it2,it2s
      real(real64) :: a2
!-----------------------------------------------
      recc = zero
      ij1=jlist(ja1)
      ij2=jlist(ja2)
      ia1=jjq1(3,ij1)-1
      ia2=jjq1(3,ij2)-1
      ib1=jjq2(3,ij1)-1
      ib2=jjq2(3,ij2)-1
      if(ja1 == 1.and.ja2 == 2) then
        it2=ia1
        it2s=ib1
        j2=jjc1(1)-1
        j2s=jjc2(1)-1
      else
        n1=ja2-1
        j2=jjc1(n1)-1
        j2s=jjc2(n1)-1
        n2=ja2-2
        it2=jjc1(n2)-1
        it2s=jjc2(n2)-1
      end if
      if(ire == 0) then
        iat = 0
        if (ittk(it2s,k1,it2) == 0) return
        if (ittk(ib2,k2,ia2) == 0) return
        if (ittk(j2s,ka,j2) == 0) return
        if (ittk(it2s,ib2,j2s) == 0) return
        if (ittk(k1,k2,ka) == 0) return
        if (ittk(it2,ia2,j2) == 0) return
        iat = 1
        return
      end if
      call nine(it2s,k1,it2,ib2,k2,ia2,j2s,ka,j2,0,iat,a2)
      recc=a2*dsqrt(dble((it2+1)*(ka+1)*(ia2+1)*(j2s+1)))
      return
      end subroutine diaga4
end module diaga4_i
