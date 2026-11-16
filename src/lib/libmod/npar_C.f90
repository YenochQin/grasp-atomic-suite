!
!***********************************************************************
!                                                                      *
      module npar_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:23:52  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: nuclear_model_parameter_count
      real(double), dimension(2) :: nuclear_model_parameters
      end module npar_c

      module nsmdat_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:36:34  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double) :: sqn,dmomnm,qmomb
      real(double) :: hfsi, hfsd, hfsq
      real(double) :: smsi, smsd, smsq
      end module nsmdat_c
