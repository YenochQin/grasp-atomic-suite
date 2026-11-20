      module initial_orbital_rotation_matrix
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer, parameter :: nlmax = 20
      real(real64) :: cici(20*nlmax*nlmax)
      end module initial_orbital_rotation_matrix
