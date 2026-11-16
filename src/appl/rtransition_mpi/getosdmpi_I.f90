      module getosd_i
      implicit none (type, external)
      interface
      subroutine getosd (name)
      implicit none (type, external)
      character (len = 128) :: name(2)
      end subroutine
      end interface
      end module
