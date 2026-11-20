      module getcid_i
      implicit none (type, external)
      interface
      subroutine getcid (isotope_file_path, radial_file_path, block_labels)
      implicit none (type, external)
      character(len=*):: isotope_file_path, radial_file_path
      character(len=8) ::  block_labels(*)
      end subroutine getcid
      end interface
      end module getcid_i
