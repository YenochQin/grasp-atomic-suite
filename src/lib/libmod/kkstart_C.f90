      module kkstart_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer, parameter :: kmax = 20
      integer :: kstart(0:kmax)
      end module kkstart_c
