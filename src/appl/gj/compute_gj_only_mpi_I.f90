      module compute_gj_only_mpi_i
      implicit none (type, external)
      interface
      subroutine compute_gj_only_mpi (gjc_diag, dgjc_diag)
      use, intrinsic :: iso_fortran_env, only: real64
      use prnt_c, only: nvec
      implicit none (type, external)
      real(real64), intent(out) :: gjc_diag(nvec), dgjc_diag(nvec)
      end subroutine compute_gj_only_mpi
      end interface
      end module compute_gj_only_mpi_i
