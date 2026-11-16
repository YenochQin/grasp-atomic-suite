      module setham_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine setham (myid, nprocs, jblock, elsto,icstrt, nelmntt,  &
                         atwinv,slf_en)
      use, intrinsic :: iso_fortran_env, only: real64, int64
      implicit none (type, external)
      integer(int64)       :: nelmntt
      integer             :: jblock, icstrt
      integer, intent(in) :: myid, nprocs
      real(real64) :: elsto, atwinv
      real(real64) :: slf_en(*)
      end subroutine setham
      end interface
      end module
