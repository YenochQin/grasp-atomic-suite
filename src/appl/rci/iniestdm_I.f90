      module iniestdm_i
      implicit none (type, external)
      interface
      subroutine iniestdm (nmax, ncf, niv, basis, hmx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in):: nmax, ncf, niv
      real(real64) :: basis(*), hmx(*)
      end subroutine  iniestdm
      end interface
      end module
