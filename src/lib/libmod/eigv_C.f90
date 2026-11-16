!
!***********************************************************************
!                                                                      *
      module eigv_c
!                                                                      *
!***********************************************************************
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  07:38:02   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(double) :: eav, eavff, eavii
      real(double), dimension(:), pointer :: eval, evalff, evalii
      real(double), dimension(:), pointer :: evec, evecff, evecii
      end module eigv_c
