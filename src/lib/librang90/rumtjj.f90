!*******************************************************************
!                                                                  *
      subroutine rumtjj(knt,jj,lq,lv,l)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 16  --------------   *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2017  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mtjj_c
      use mtjj2_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use jthn_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: jj, knt
      integer, intent(out) :: lq, lv, l
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: kt, kntmin
!-----------------------------------------------

      if(jj < 9) then
        kt=mt(knt)
      else if(jj == 9) then
        if(knt > 300) then
          kntmin=knt-300
          kt=mt9(kntmin)
        else
          print*, "ERROR in RUMTJJ"
          stop
!GG          CALL RUMT67(KNT,NR,LQ,LS,L)
!GG          RETURN
        endif
      else
        kt=mt11(knt)
      endif
      lq=jthn(kt,3,100)
      lv=jthn(kt,2,100)
      l=jthn(kt,1,100)
      return
      end subroutine rumtjj
