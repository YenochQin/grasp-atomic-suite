      module kkstart_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: kmax = 20
      integer :: kstart(0:kmax)
      end module kkstart_c
