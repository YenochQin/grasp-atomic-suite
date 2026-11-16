      module rmcdhf_mem_c
!   Written by  G. Gaigalas               Vilnius, September 2021  *
      use vast_kind_param, only:  double, long
      real(double), dimension(:), pointer :: coeff_31
      real(double), dimension(:), pointer :: coeff_32
      real(double), dimension(:), pointer :: coeff_33
      real(double), dimension(:), pointer :: coeff_34
      real(double), dimension(:), pointer :: coeff_35
      real(double), dimension(:), pointer :: coeff_36
      real(double), dimension(:), pointer :: coeff_37
      real(double), dimension(:), pointer :: coeff_38
      real(double), dimension(:), pointer :: coeff_39
      real(double), dimension(:), pointer :: coeff_40
      real(double), dimension(:), pointer :: coeff_41
      real(double), dimension(:), pointer :: coeff_42
      real(double), dimension(:), pointer :: coeff_43
      real(double), dimension(:), pointer :: coeff_44
      real(double), dimension(:), pointer :: coeff_45
      real(double), dimension(:), pointer :: coeff_46
      real(double), dimension(:), pointer :: coeff_47
      real(double), dimension(:), pointer :: coeff_48
      real(double), dimension(:), pointer :: coeff_49
      real(double), dimension(:), pointer :: coeff_50
      real(double), dimension(:), pointer :: coeff
      integer, dimension(:), pointer ::imin1_30, imax1_30
      integer, dimension(:), pointer ::imin2_30, imax2_30
      integer, dimension(:), pointer ::indx, iclmn
      integer, dimension(:), pointer ::iendc_30, irow_30
      integer, dimension(:), pointer ::lab_31,ncontr_31,iclmn_31,indx_31
      integer, dimension(:), pointer ::lab_32,ncontr_32,iclmn_32,indx_32
      integer, dimension(:), pointer ::lab_33,ncontr_33,iclmn_33,indx_33
      integer, dimension(:), pointer ::lab_34,ncontr_34,iclmn_34,indx_34
      integer, dimension(:), pointer ::lab_35,ncontr_35,iclmn_35,indx_35
      integer, dimension(:), pointer ::lab_36,ncontr_36,iclmn_36,indx_36
      integer, dimension(:), pointer ::lab_37,ncontr_37,iclmn_37,indx_37
      integer, dimension(:), pointer ::lab_38,ncontr_38,iclmn_38,indx_38
      integer, dimension(:), pointer ::lab_39,ncontr_39,iclmn_39,indx_39
      integer, dimension(:), pointer ::lab_40,ncontr_40,iclmn_40,indx_40
      integer, dimension(:), pointer ::lab_41,ncontr_41,iclmn_41,indx_41
      integer, dimension(:), pointer ::lab_42,ncontr_42,iclmn_42,indx_42
      integer, dimension(:), pointer ::lab_43,ncontr_43,iclmn_43,indx_43
      integer, dimension(:), pointer ::lab_44,ncontr_44,iclmn_44,indx_44
      integer, dimension(:), pointer ::lab_45,ncontr_45,iclmn_45,indx_45
      integer, dimension(:), pointer ::lab_46,ncontr_46,iclmn_46,indx_46
      integer, dimension(:), pointer ::lab_47,ncontr_47,iclmn_47,indx_47
      integer, dimension(:), pointer ::lab_48,ncontr_48,iclmn_48,indx_48
      integer, dimension(:), pointer ::lab_49,ncontr_49,iclmn_49,indx_49
      integer, dimension(:), pointer ::lab_50,ncontr_50,iclmn_50,indx_50
      end module rmcdhf_mem_c
