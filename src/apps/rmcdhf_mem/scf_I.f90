      module scf_i
      implicit none (type, external)
      interface
      subroutine scf (extended_optimal_level, output_wavefunction_path)
      implicit none (type, external)
      logical, intent(in) :: extended_optimal_level
      character (len = *) :: output_wavefunction_path
      end subroutine scf
      end interface
      end module scf_i
