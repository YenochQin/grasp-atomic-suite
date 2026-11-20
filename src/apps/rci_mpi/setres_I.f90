      module setres_i
      implicit none (type, external)
      interface
      subroutine setres (isotope_file_path, radial_file_path, block_labels)
      implicit none (type, external)
      character (len = *) :: isotope_file_path
      character (len = *) :: radial_file_path
      character (len = 8) :: block_labels(*)
      end subroutine setres
      end interface
      end module setres_i
