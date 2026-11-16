      module genintrk_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine genintrk (myid, nprocs, n, j2max)
        implicit none (type, external)
        integer, intent(in) :: myid, nprocs
        integer, intent(out) :: n, j2max
      end subroutine
      end interface
      end module
