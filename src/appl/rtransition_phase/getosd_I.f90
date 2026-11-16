      module getosd_i
      implicit none (type, external)
      interface
      subroutine getosd (isofile,name)
      implicit none (type, external)
      character(len=*):: isofile
      character (len = 24) :: name(2)
      end subroutine getosd
      end interface
      end module getosd_i
