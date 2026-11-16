      module setham_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine setham (myid, nprocs, jblock, elsto,icstrt, nelmntt,  &
                         atwinv,slf_en)
      use vast_kind_param,   only: double, long
      integer(long)       :: nelmntt
      integer             :: jblock, icstrt
      integer, intent(in) :: myid, nprocs
      real(double) :: elsto, atwinv
      real(double), dimension(*) :: slf_en
      end subroutine setham
      end interface
      end module
