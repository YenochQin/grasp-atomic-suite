!*******************************************************************
!                                                                  *
      integer function itjj(ik,id,kg,bk,bd,ibt,bt,kg1,itp,itg,iq)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 06  --------------   *
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
      integer,      intent(in)                :: kg, iq
      integer,      intent(out)               :: kg1, itp, itg
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
      itjj=0
      if(id(3) > 37) return
      kg1=2*kg
      if(ittk(id(6),ik(6),kg1) == 0)return
      itk=ik(1)
      itd=id(1)
      if(id(3) < 9) then
        itp1=imptjj(itk)
        itp=imptjj(itd)
        if(itp1 /= itp)return
        itg1=imgtjj(itk)
        itg=imgtjj(itd)
      elseif(id(3) == 9) then
        if(itk > 300) then
          if(itd < 300) call mes(51)
          if(id(4) > 2) call mes(11)
          if(ik(4) > 2) call mes(11)
          itk=itk-300
          itd=itd-300
          itp1=imptjj9(itk)
          itp=imptjj9(itd)
          if(itp1 /= itp)return
          itg1=imgtjj9(itk)
          itg=imgtjj9(itd)
        else
          print*, "ERROR in ITJJ"
          stop
        end if
      else
        if(id(4) > 2) call mes(11)
        if(ik(4) > 2) call mes(11)
        itp1=imptjj11(itk)
        itp=imptjj11(itd)
        if(itp1 /= itp)return
        itg1=imgtjj11(itk)
        itg=imgtjj11(itd)
      endif
      if(itg1 /= itg)return
      itjj=1
      ibt(2)=id(2)
      ibt(3)=id(3)
      ibt(4)=id(4)+iq
      bt(3)=bd(3)+half*dble(iq)
      return
      end function itjj
