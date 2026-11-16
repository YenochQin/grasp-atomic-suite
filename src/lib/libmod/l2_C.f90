      module l2_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer, parameter :: mtriad = 12
      integer :: j2s(mtriad, 3), j3s(mtriad, 3)
      end module l2_c
