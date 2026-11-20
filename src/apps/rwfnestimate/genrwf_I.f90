      module genrwf_i
      implicit none (type, external)
      interface
      subroutine genrwf(default_orbitals)
      implicit none (type, external)
      character(len=*), intent(in) :: default_orbitals
      end subroutine genrwf
      end interface
      end module genrwf_i
