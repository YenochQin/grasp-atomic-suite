      module alloccheck_i
      implicit none (type, external)
      interface
      subroutine alloccheck(n, irez)
      implicit none (type, external)
      integer, intent(in)  :: n
      integer, intent(out) :: irez
      end subroutine alloccheck
      end interface
      end module alloccheck_i
