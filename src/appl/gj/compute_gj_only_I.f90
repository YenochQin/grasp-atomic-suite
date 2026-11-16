      module compute_gj_only_i
      interface
      subroutine compute_gj_only (gjc_diag, dgjc_diag)
      use vast_kind_param, only: double
      use prnt_c, only: nvec
      real(double), intent(out) :: gjc_diag(nvec), dgjc_diag(nvec)
      end subroutine
      end interface
      end module
