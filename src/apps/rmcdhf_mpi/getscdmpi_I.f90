      module getscdmpi_i
      implicit none (type, external)
      interface
      subroutine getscdmpi (extended_optimal_level, block_labels, isotope_file_path, radial_file_path)
      implicit none (type, external)
      logical, intent(out) :: extended_optimal_level
      character (len = 8) :: block_labels(*)
      character (len = *) :: isotope_file_path
      character (len = *) :: radial_file_path
      end subroutine getscdmpi
      end interface
      end module getscdmpi_i
