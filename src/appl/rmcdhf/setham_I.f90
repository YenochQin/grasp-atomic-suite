      module setham_i
      implicit none (type, external)
      interface
      subroutine setham (jblock, myid, nprocs)
      implicit none (type, external)
      integer, intent(in) :: jblock
      integer, intent(in) :: myid
      integer, intent(in) :: nprocs
      end subroutine
      end interface
      end module
