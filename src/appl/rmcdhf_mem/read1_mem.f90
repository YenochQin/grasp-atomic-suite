!*******************************************************************
!                                                                  *
      subroutine read1_mem(nfile, ncontr_tot2, lab, ncontr)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: int64
      use rmcdhf_mem_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nfile
      integer :: lab, ncontr
      integer(int64), intent(in) :: ncontr_tot2
!-----------------------------------------------

      if(nfile == 31) then
         lab    = lab_31(ncontr_tot2)
         ncontr = ncontr_31(ncontr_tot2)
      else if(nfile == 32) then
         lab    = lab_32(ncontr_tot2)
         ncontr = ncontr_32(ncontr_tot2)
      else if(nfile == 33) then
         lab    = lab_33(ncontr_tot2)
         ncontr = ncontr_33(ncontr_tot2)
      else if(nfile == 34) then
         lab    = lab_34(ncontr_tot2)
         ncontr = ncontr_34(ncontr_tot2)
      else if(nfile == 35) then
         lab    = lab_35(ncontr_tot2)
         ncontr = ncontr_35(ncontr_tot2)
      else if(nfile == 36) then
         lab    = lab_36(ncontr_tot2)
         ncontr = ncontr_36(ncontr_tot2)
      else if(nfile == 37) then
         lab    = lab_37(ncontr_tot2)
         ncontr = ncontr_37(ncontr_tot2)
      else if(nfile == 38) then
         lab    = lab_38(ncontr_tot2)
         ncontr = ncontr_38(ncontr_tot2)
      else if(nfile == 39) then
         lab    = lab_39(ncontr_tot2)
         ncontr = ncontr_39(ncontr_tot2)
      else if(nfile == 40) then
         lab    = lab_40(ncontr_tot2)
         ncontr = ncontr_40(ncontr_tot2)
      else if(nfile == 41) then
         lab    = lab_41(ncontr_tot2)
         ncontr = ncontr_41(ncontr_tot2)
      else if(nfile == 42) then
         lab    = lab_42(ncontr_tot2)
         ncontr = ncontr_42(ncontr_tot2)
      else if(nfile == 43) then
         lab    = lab_43(ncontr_tot2)
         ncontr = ncontr_43(ncontr_tot2)
      else if(nfile == 44) then
         lab    = lab_44(ncontr_tot2)
         ncontr = ncontr_44(ncontr_tot2)
      else if(nfile == 45) then
         lab    = lab_45(ncontr_tot2)
         ncontr = ncontr_45(ncontr_tot2)
      else if(nfile == 46) then
         lab    = lab_46(ncontr_tot2)
         ncontr = ncontr_46(ncontr_tot2)
      else if(nfile == 47) then
         lab    = lab_47(ncontr_tot2)
         ncontr = ncontr_47(ncontr_tot2)
      else if(nfile == 48) then
         lab    = lab_48(ncontr_tot2)
         ncontr = ncontr_48(ncontr_tot2)
      else if(nfile == 49) then
         lab    = lab_49(ncontr_tot2)
         ncontr = ncontr_49(ncontr_tot2)
      else if(nfile == 50) then
         lab    = lab_50(ncontr_tot2)
         ncontr = ncontr_50(ncontr_tot2)
      else
         print*, "Error in read1_memi  NFILE=",nfile
         stop
      end if

      return
      end subroutine read1_mem
