      module orbout_i
      implicit none (type, external)
      interface
      subroutine orbout (output_wavefunction_path)
      implicit none (type, external)
      character (len = *), intent(in) :: output_wavefunction_path
      end subroutine orbout
      end interface
      end module orbout_i
