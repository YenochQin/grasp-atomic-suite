      module fco_i
      implicit none (type, external)
      interface
      real(real64) function fco (multipole_rank, csf_index, orbital_a, orbital_b)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: multipole_rank
      integer :: csf_index
      integer, intent(in) :: orbital_a
      integer, intent(in) :: orbital_b
      end function fco
      end interface
      end module fco_i
