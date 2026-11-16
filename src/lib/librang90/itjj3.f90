module itjj3_i
    implicit none (type, external)
    private
    public :: itjj3
contains
!******************************************************************
!                                                                  *
      integer function itjj3(ik,id,kg1,bk,bd,ibt,bt,itp,itg,iq)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 08  --------------   *
!                                                                  *
!     FUNCTION CALLED: ITTK                                        *
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
      use cons_c,          only: half
      use ribojj_c
      use ribojj9_c
      use ribojj11_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use ittk_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)                :: kg1, iq
      integer,      intent(out)               :: itp, itg
      integer,      intent(in) :: ik(7), id(7)
      integer,      intent(out) :: ibt(7)
      real(real64), intent(in) :: bk(3), bd(3)
      real(real64), intent(out) :: bt(3)
!      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3),BK(3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: itk, itd, itp1, itg1
!-----------------------------------------------
      itjj3=0
      if(id(3) > 37) return
      if(ittk(id(6),ik(6),kg1) == 0)return
      itk=ik(1)
      itd=id(1)
      if(id(3) < 9) then
        itp1=imptjj(itk)
        itp=impnjj(itd)
        if(itp1 /= itp)return
        itg1=imgtjj(itk)
        itg=imgnjj(itd)
      elseif(id(3) == 9) then
        if(itk > 300) then
          if(itd < 300) call mes(53)
          if(id(4) > 2) call mes(13)
          if(ik(4) > 2) call mes(13)
          itk=itk-300
          itd=itd-300
          itp1=imptjj9(itk)
          itp=impnjj9(itd)
          if(itp1 /= itp)return
          itg1=imgtjj9(itk)
          itg=imgnjj9(itd)
        else
          print*, "ERROR in ITJJ3"
          stop
        endif
      else
        if(id(4) > 2) call mes(13)
        if(ik(4) > 2) call mes(13)
        itp1=imptjj11(itk)
        itp=impnjj11(itd)
        if(itp1 /= itp)return
        itg1=imgtjj11(itk)
        itg=imgnjj11(itd)
      endif
      if(itg1 /= itg)return
      itjj3=1
      ibt(2)=id(2)
      ibt(3)=id(3)
      ibt(4)=id(4)+iq
      bt(3)=bd(3)+half*dble(iq)
      return
      end function itjj3
end module itjj3_i
