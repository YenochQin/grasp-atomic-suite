!*******************************************************************
!                                                                  *
      integer function itjj2(ik,id,kg1,bk,bd,ibt,bt,itp,itg,iq)
!                                                                  *
!   ---------------  SECTION SQJJ  SUBPROGRAM 07  --------------   *
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
      use vast_kind_param, only: double
      use cons_c,          only: half
      use ribojj_c
      use ribojj9_c
      use ribojj11_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use mes_i
      use ittk_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer,      intent(in)                :: kg1, iq
      integer,      intent(out)               :: itp, itg
      integer,      intent(in),  dimension(7) :: ik, id
      integer,      intent(out), dimension(7) :: ibt
      real(double), intent(in),  dimension(3) :: bk, bd
      real(double), intent(out), dimension(3) :: bt
!      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3),BK(3)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: itk, itd, itp1, itg1
!-----------------------------------------------
      itjj2=0
      if(id(3) > 37) return
      if(ittk(id(6),ik(6),kg1) == 0)return
      itk=ik(1)
      itd=id(1)
      if(id(3) < 9) then
        itp1=impnjj(itk)
        itp=imptjj(itd)
        if(itp1.ne.itp)return
        itg1=imgnjj(itk)
        itg=imgtjj(itd)
      elseif(id(3) == 9) then
        if(itk > 300) then
          if(itd < 300) call mes(52)
          if(id(4) > 2) call mes(12)
          if(ik(4) > 2) call mes(12)
          itk=itk-300
          itd=itd-300
          itp1=impnjj9(itk)
          itp=imptjj9(itd)
          if(itp1 /= itp)return
          itg1=imgnjj9(itk)
          itg=imgtjj9(itd)
        else
          print*, "ERROR in ITJJ2"
          stop
        end if
      else
        if(id(4) > 2) call mes(12)
        if(ik(4) > 2) call mes(12)
        itp1=impnjj11(itk)
        itp=imptjj11(itd)
        if(itp1 /= itp)return
        itg1=imgnjj11(itk)
        itg=imgtjj11(itd)
      endif
      if(itg1 /= itg)return
      itjj2=1
      ibt(2)=id(2)
      ibt(3)=id(3)
      ibt(4)=id(4)+iq
      bt(3)=bd(3)+half*dble(iq)
      return
      end function itjj2
