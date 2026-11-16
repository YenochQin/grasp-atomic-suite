      module setisompi_i
      implicit none (type, external)
      interface
      subroutine setisompi (isofile)
      implicit none (type, external)
      character(len=*) :: isofile
      end subroutine
      end interface
      end module
