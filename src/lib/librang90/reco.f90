module reco_i
    implicit none (type, external)
    private
    public :: reco
contains
!*******************************************************************
!                                                                  *
      subroutine reco(ja1,ja2,ja3,ja4,ka,iat)
!                                                                  *
!     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
!                                                                  *
!     NO SUBROUTINE CALLED                                         *
!                                                                  *
!   Written by  G. Gaigalas                                        *
!   Transform to fortran 90/95 by G. Gaigalas       December 2012  *
!   The last modification made by G. Gaigalas       October  2020  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use m_c,             only: npeel, jlist, jjc1, jjc2, jjq1, jjq2
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: ja1, ja2, ja3, ja4, ka
      integer, intent(out) :: iat
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ij1, ij2, ij, j, ia1, ia2
!-----------------------------------------------
      iat=1
      if(npeel == 1)return
      ij1=jlist(ja1)
      ij2=jlist(ja2)
      if(ja1 == 1.and.ja2 == 2) go to 1
      if(ka /= 0)go to 5
!
!  CASES WHEN :          KA = 0
!                  OR    JA1 = JA2
!                  OR    JA1 = 1    JA2 = 2
!
    1 do i=1,npeel
        ij=jlist(i)
        if(i < npeel-1) then
          if(jjc1(i) /= jjc2(i))iat=0
        end if
        if(ka /= 0) then
          if(i == ja1) cycle
          if(i == ja2) cycle
        end if
        do j=1,3
          if(ja1 == ja2 .and. i == ja1 .and. j == 1) cycle
          if(jjq1(j,ij) /= jjq2(j,ij))iat=0
        end do
      end do
      return
!
!  OTHER CASES
!
    5 continue
      do i=1,npeel
        ij=jlist(i)
        if(i < npeel-1) then
          ia1=ja1-1
          ia2=ja2-1
          if(ja1 == 1)ia1=ja1
          if(i >= ia1.and.i < ia2)go to 7
          if(jjc1(i) /= jjc2(i))iat=0
        end if
    7   if(i == ja1) cycle
        if(i == ja2) cycle
        if((ka == 2).and.(i == ja3)) cycle
        if((ka == 3).and.(i == ja3)) cycle
        if((ka == 3).and.(i == ja4)) cycle
        do j=1,3
          if(jjq1(j,ij) /= jjq2(j,ij))iat=0
        end do
      end do
      return
      end subroutine reco
end module reco_i
