      module cimat_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
      integer, parameter :: nlmax = 20
      real(real64) :: cfci(20*nlmax*nlmax)
      real(real64) :: cici(20*nlmax*nlmax)
      end module cimat_c
