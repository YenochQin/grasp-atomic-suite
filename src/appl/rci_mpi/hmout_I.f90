      module hmout_i
      implicit none (type, external)
      interface
      subroutine hmout (imcdf)
      implicit none (type, external)
      integer, intent(in) :: imcdf
      end subroutine
      end interface
      end module
