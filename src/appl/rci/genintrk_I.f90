      module genintrk_i
      implicit none (type, external)
      interface
      subroutine genintrk (myid, nprocs, n, j2max)
        implicit none (type, external)
        integer, intent(in) :: myid, nprocs
        integer, intent(out) :: n, j2max
      end subroutine
      end interface
      end module
