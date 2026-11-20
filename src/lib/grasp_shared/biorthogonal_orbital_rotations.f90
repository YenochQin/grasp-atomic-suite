      module biorthogonal_orbital_rotations
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer, parameter :: nlmax = 20
      real(real64) :: cfci(20*nlmax*nlmax)
      real(real64) :: cici(20*nlmax*nlmax)
      end module biorthogonal_orbital_rotations
