      module hmout_i
      implicit none (type, external)
      interface
      subroutine hmout (imcdf)
      implicit none (type, external)
      integer, intent(in) :: imcdf
      end subroutine hmout
      end interface
      end module hmout_i
