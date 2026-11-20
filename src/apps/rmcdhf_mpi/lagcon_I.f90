      module lagcon_i
      implicit none (type, external)
      interface
      subroutine lagcon (orbital_index, mpi_process_count)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer, intent(in) :: mpi_process_count
      end subroutine lagcon
      end interface
      end module lagcon_i
