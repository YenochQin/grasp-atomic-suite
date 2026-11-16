      module couple_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
      integer, parameter :: mangm = 60
      integer, parameter :: mtriad = 12
      integer :: j1(mangm)
      integer :: j2(mtriad, 3), j3(mtriad, 3)
      integer :: mja, nja
      logical :: free(mangm)
      end module couple_c
