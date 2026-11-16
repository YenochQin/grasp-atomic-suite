      module kapdata_i
      implicit none (type, external)
      interface
      subroutine kapdata (ntestg, ncore1, ncore2)
      implicit none (type, external)
      integer, intent(in) :: ntestg
      integer, intent(in) :: ncore1
      integer, intent(in) :: ncore2
      end subroutine
      end interface
      end module
