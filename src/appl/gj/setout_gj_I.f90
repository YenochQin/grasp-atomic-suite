      module setout_gj_i
      implicit none (type, external)
      interface
      subroutine setout_gj (name, nci, gjc_diag, dgjc_diag)
      use, intrinsic :: iso_fortran_env, only: real64
      use prnt_c, only: nvec
      implicit none (type, external)
      character(len=24), intent(in) :: name
      integer, intent(in) :: nci
      real(real64), intent(in) :: gjc_diag(nvec), dgjc_diag(nvec)
      end subroutine setout_gj
      end interface
      end module setout_gj_i
