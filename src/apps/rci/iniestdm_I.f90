      module iniestdm_i
      implicit none (type, external)
      interface
      subroutine iniestdm (nmax, csf_count, initial_vector_count, basis, hmx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in):: nmax, csf_count, initial_vector_count
      real(real64) :: basis(*), hmx(*)
      end subroutine  iniestdm
      end interface
      end module iniestdm_i
