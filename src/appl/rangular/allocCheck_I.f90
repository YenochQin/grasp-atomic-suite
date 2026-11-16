      module alloccheck_i
      implicit none (type, external)
      interface
      subroutine alloccheck(n, irez)
      implicit none (type, external)
      integer, intent(in)  :: n
      integer, intent(out) :: irez
      end subroutine
      end interface
      end module
