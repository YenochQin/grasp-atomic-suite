      module scfmpi_i
      implicit none (type, external)
      interface
      subroutine scfmpi (extended_optimal_level, output_wavefunction_path)
      implicit none (type, external)
      logical, intent(in) :: extended_optimal_level
      character (len = *) :: output_wavefunction_path
      end subroutine scfmpi
      end interface
      end module scfmpi_i
