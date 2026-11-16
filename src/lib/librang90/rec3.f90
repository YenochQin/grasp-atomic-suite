module rec3_i
    implicit none (type, external)
    private
    public :: rec3
contains
!*******************************************************************
!                                                                  *
      subroutine rec3(ja1,ja2,ja3,k1,k2,k3,ire,iat,recc)
!                                                                  *
!   ---------------  SECTION REC    SUBPROGRAM 07  --------------  *
!                                                                  *
!     SUBROUTINE CALLED:  RECO3                                    *
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
      use, intrinsic :: iso_fortran_env, only: real64
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use reco3_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)       :: ja1,ja2,ja3,k1,k2,k3,ire
      integer, intent(out)      :: iat
      real(real64), intent(out) :: recc
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ifaz
!-----------------------------------------------
      if((ja3 > ja1).and.(ja3 > ja2)) then
        if(ja1-ja2 < 0) then
          call reco3(ja1,ja2,ja3,k1,k2,k3,ire,iat,recc)
        else if(ja1-ja2 > 0) then
          call reco3(ja2,ja1,ja3,k2,k1,k3,ire,iat,recc)
          ifaz=k1+k2-k3
          if((ifaz/4)*4 /= ifaz)recc=-recc
        else
          go to 10
        end if
      else if((ja3 < ja1).and.(ja3 < ja2)) then
        if(ja1-ja2 < 0) then
          call reco3(ja3,ja1,ja2,k3,k1,k2,ire,iat,recc)
          if((k3/2)*2 /= k3)recc=-recc
        else if(ja1-ja2 > 0) then
          call reco3(ja3,ja2,ja1,k3,k2,k1,ire,iat,recc)
          ifaz=k1+k2+k3
          if((ifaz/4)*4 /= ifaz)recc=-recc
        else
          go to 10
        end if
      else
        if(ja1-ja2 < 0)then
          call reco3(ja1,ja3,ja2,k1,k3,k2,ire,iat,recc)
          ifaz=k1-k2-k3
          if((ifaz/4)*4 /= ifaz)recc=-recc
        else if(ja1-ja2 > 0) then
          call reco3(ja2,ja3,ja1,k2,k3,k1,ire,iat,recc)
          if((k1/2)*2 /= k1)recc=-recc
        else
          go to 10
        end if
      end if
      return
   10 write(99,100)
  100 format(5x,'ERRO IN REC')
      stop
      end subroutine rec3
end module rec3_i
