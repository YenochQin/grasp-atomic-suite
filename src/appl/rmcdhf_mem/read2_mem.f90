!*******************************************************************
!                                                                  *
      subroutine read2_mem(nfile,ncontr_tot,ncontr,iclmndum)
!                                                                  *
!   Written by  G. Gaigalas               Vilnius, September 2021  *
!                                                                  *
!*******************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: int64
      use rmcdhf_mem_c, only: coeff, coeff_31, coeff_32, coeff_33, coeff_34, coeff_35, coeff_36, coeff_37, coeff_38, &
         & coeff_39, coeff_40, coeff_41, coeff_42, coeff_43, coeff_44, coeff_45, coeff_46, coeff_47, coeff_48, &
         & coeff_49, coeff_50, iclmn_31, iclmn_32, iclmn_33, iclmn_34, iclmn_35, iclmn_36, iclmn_37, iclmn_38, &
         & iclmn_39, iclmn_40, iclmn_41, iclmn_42, iclmn_43, iclmn_44, iclmn_45, iclmn_46, iclmn_47, iclmn_48, &
         & iclmn_49, iclmn_50, indx, indx_31, indx_32, indx_33, indx_34, indx_35, indx_36, indx_37, indx_38, &
         & indx_39, indx_40, indx_41, indx_42, indx_43, indx_44, indx_45, indx_46, indx_47, indx_48, indx_49, indx_50
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer(int64), intent(in) :: ncontr_tot
      integer, intent(in) :: nfile
      integer :: ncontr, iclmndum
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------

      do i = 1, ncontr
      !*** Read column index, sparse matrix index, and
      !*** coefficient for all contributions from this integral
         if(nfile == 31) then
            iclmndum = iclmn_31(ncontr_tot+i)
            indx(i) = indx_31(ncontr_tot+i)
            coeff(i) = coeff_31(ncontr_tot+i)
         else if(nfile == 32) then
            iclmndum = iclmn_32(ncontr_tot+i)
            indx(i) = indx_32(ncontr_tot+i)
            coeff(i) = coeff_32(ncontr_tot+i)
         else if(nfile == 33) then
            iclmndum = iclmn_33(ncontr_tot+i)
            indx(i) = indx_33(ncontr_tot+i)
            coeff(i) = coeff_33(ncontr_tot+i)
         else if(nfile == 34) then
            iclmndum = iclmn_34(ncontr_tot+i)
            indx(i) = indx_34(ncontr_tot+i)
            coeff(i) = coeff_34(ncontr_tot+i)
         else if(nfile == 35) then
            iclmndum = iclmn_35(ncontr_tot+i)
            indx(i) = indx_35(ncontr_tot+i)
            coeff(i) = coeff_35(ncontr_tot+i)
         else if(nfile == 36) then
            iclmndum = iclmn_36(ncontr_tot+i)
            indx(i) = indx_36(ncontr_tot+i)
            coeff(i) = coeff_36(ncontr_tot+i)
         else if(nfile == 37) then
            iclmndum = iclmn_37(ncontr_tot+i)
            indx(i) = indx_37(ncontr_tot+i)
            coeff(i) = coeff_37(ncontr_tot+i)
         else if(nfile == 38) then
            iclmndum = iclmn_38(ncontr_tot+i)
            indx(i) = indx_38(ncontr_tot+i)
            coeff(i) = coeff_38(ncontr_tot+i)
         else if(nfile == 39) then
            iclmndum = iclmn_39(ncontr_tot+i)
            indx(i) = indx_39(ncontr_tot+i)
            coeff(i) = coeff_39(ncontr_tot+i)
         else if(nfile == 40) then
            iclmndum = iclmn_40(ncontr_tot+i)
            indx(i) = indx_40(ncontr_tot+i)
            coeff(i) = coeff_40(ncontr_tot+i)
         else if(nfile == 41) then
            iclmndum = iclmn_41(ncontr_tot+i)
            indx(i) = indx_41(ncontr_tot+i)
            coeff(i) = coeff_41(ncontr_tot+i)
         else if(nfile == 42) then
            iclmndum = iclmn_42(ncontr_tot+i)
            indx(i) = indx_42(ncontr_tot+i)
            coeff(i) = coeff_42(ncontr_tot+i)
         else if(nfile == 43) then
            iclmndum = iclmn_43(ncontr_tot+i)
            indx(i) = indx_43(ncontr_tot+i)
            coeff(i) = coeff_43(ncontr_tot+i)
         else if(nfile == 44) then
            iclmndum = iclmn_44(ncontr_tot+i)
            indx(i) = indx_44(ncontr_tot+i)
            coeff(i) = coeff_44(ncontr_tot+i)
         else if(nfile == 45) then
            iclmndum = iclmn_45(ncontr_tot+i)
            indx(i) = indx_45(ncontr_tot+i)
            coeff(i) = coeff_45(ncontr_tot+i)
         else if(nfile == 46) then
            iclmndum = iclmn_46(ncontr_tot+i)
            indx(i) = indx_46(ncontr_tot+i)
            coeff(i) = coeff_46(ncontr_tot+i)
         else if(nfile == 47) then
            iclmndum = iclmn_47(ncontr_tot+i)
            indx(i) = indx_47(ncontr_tot+i)
            coeff(i) = coeff_47(ncontr_tot+i)
         else if(nfile == 48) then
            iclmndum = iclmn_48(ncontr_tot+i)
            indx(i) = indx_48(ncontr_tot+i)
            coeff(i) = coeff_48(ncontr_tot+i)
         else if(nfile == 49) then
            iclmndum = iclmn_49(ncontr_tot+i)
            indx(i) = indx_49(ncontr_tot+i)
            coeff(i) = coeff_49(ncontr_tot+i)
         else if(nfile == 50) then
            iclmndum = iclmn_50(ncontr_tot+i)
            indx(i) = indx_50(ncontr_tot+i)
            coeff(i) = coeff_50(ncontr_tot+i)
         else
            print*, "Error in read2_mem"
            stop
         end if
      end do

      return
      end subroutine read2_mem
