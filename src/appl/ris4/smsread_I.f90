      module smsread_i
      interface
!...Translated by Gediminas Gaigalas 11/18/19
      subroutine smsread (vint,vint2)
      use vast_kind_param, only: double
      use parameter_def,   only: nnnw
      real(double), dimension(nnnw,nnnw), intent(in) :: vint, vint2
      end subroutine
      end interface
      end module
