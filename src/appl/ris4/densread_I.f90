      module densread_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine densread (dint1,dint2,dint3,                          &
                           dint4,dint5,dint6,                          &
                           dint7)
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw, nnnp
      use prnt_c,           only: nvec
      real(double), dimension(nnnw,nnnw), intent(in) :: dint1, dint2,  &
                                          dint3, dint4, dint5, dint6,  &
                                          dint7
      end subroutine
      end interface
      end module
