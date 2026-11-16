      module setout_gj_i
      interface
      subroutine setout_gj (name, nci, gjc_diag, dgjc_diag)
      use vast_kind_param, only: double
      use prnt_c, only: nvec
      character(len=24), intent(in) :: name
      integer, intent(in) :: nci
      real(double), intent(in) :: gjc_diag(nvec), dgjc_diag(nvec)
      end subroutine
      end interface
      end module
