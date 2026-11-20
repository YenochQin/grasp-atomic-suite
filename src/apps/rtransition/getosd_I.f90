      module getosd_i
      implicit none (type, external)
      interface
      subroutine getosd (isotope_file_path,name)
      implicit none (type, external)
      character(len=*):: isotope_file_path
      character (len = 24) :: name(2)
      end subroutine getosd
      end interface
      end module getosd_i
