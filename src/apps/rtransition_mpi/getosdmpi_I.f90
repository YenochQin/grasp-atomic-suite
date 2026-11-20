      module getosd_i
      implicit none (type, external)
      interface
      subroutine getosd (name)
      implicit none (type, external)
      character (len = 128) :: name(2)
      end subroutine getosd
      end interface
      end module getosd_i
