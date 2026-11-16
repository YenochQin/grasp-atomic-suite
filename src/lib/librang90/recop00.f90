module recop00_i
    implicit none (type, external)
    private
    public :: recop00
contains
!*******************************************************************
!                                                                  *
      subroutine recop00(ns,ja1,ja2,ka,iat)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 06  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use m_c,             only: npeel, jjc1, jjc2, jlist, jjq1, jjq2
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use ittk_i, only: ittk
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ns, ja1, ja2, ka
      integer, intent(out) :: iat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, iskr, ij1, j, ji, jj, npeelgg, kk, lk, ld
!-----------------------------------------------
      iat=1
      if(npeel == 1 .and. ns == -1)return
      iat=0
      if(ns == -1) then
         npeelgg = npeel
      else
         npeelgg = ns
      end if
      if(npeelgg > 1) then
        lk=jjc1(npeelgg-1)-1
        ld=jjc2(npeelgg-1)-1
      else if(npeelgg == 1) then
        lk=jjc1(npeelgg)-1
        ld=jjc2(npeelgg)-1
      else
        print*, "ERROR in RECOP00"
        stop
      end if
      if(ittk(lk,ld,2*ka) == 0)return
      iat=1
      if(npeel == 1)return
      do i=1,npeel
        if(ja1 /= i) then
          if(ja2 /= i) then
            ij1=jlist(i)
            if(jjq1(1,ij1) /= jjq2(1,ij1))iat=0
            if(jjq1(2,ij1) /= jjq2(2,ij1))iat=0
            if(jjq1(3,ij1) /= jjq2(3,ij1))iat=0
          end if
        end if
      end do
      if(iat == 0)return
      if(npeelgg <= 2)return
      if(ja1 <= 2)return
      do j=3,ja1
        jj=j-2
        if(jjc1(jj) /= jjc2(jj))iat=0
        if(iat == 0)return
      end do
      iskr=npeelgg-ja2
      if(iskr > 0) then
        do ji=1,iskr
          kk=ja2-2+ji
          lk=jjc1(kk)-1
          ld=jjc2(kk)-1
          if(ittk(lk,ld,2*ka) == 0)iat=0
          if(iat == 0)return
        end do
      end if
      return
      end subroutine recop00
end module recop00_i
