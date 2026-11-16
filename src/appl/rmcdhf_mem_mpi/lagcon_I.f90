      module lagcon_i
      implicit none (type, external)
      interface
      subroutine lagcon (j, nprocs)
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(in) :: nprocs
      end subroutine
      end interface
      end module
