      module cimat_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: nlmax = 20
      real(real64) :: cfci(20*nlmax*nlmax)
      real(real64) :: cici(20*nlmax*nlmax)
      end module cimat_c
