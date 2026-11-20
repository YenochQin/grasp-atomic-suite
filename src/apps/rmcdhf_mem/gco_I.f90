      module gco_i
      implicit none (type, external)
      interface
      real(real64) function gco (multipole_rank, csf_index, orbital_a, orbital_b)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer :: multipole_rank
      integer :: csf_index
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_b
      end function gco
      end interface
      end module gco_i
