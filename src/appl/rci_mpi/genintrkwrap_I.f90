      module genintrkwrap_i
      implicit none (type, external)
      interface
      subroutine genintrkwrap (myid, nprocs, j2max)
      implicit none (type, external)
      integer :: myid
      integer :: nprocs
      integer :: j2max
      end subroutine
      end interface
      end module
