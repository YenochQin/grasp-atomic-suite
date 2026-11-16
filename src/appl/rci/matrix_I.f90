      module matrix_i
      implicit none (type, external)
      interface
      subroutine matrix (ncore, j2max)
      implicit none (type, external)
      integer,  intent(in):: ncore, j2max
      end subroutine
      end interface
      end module
