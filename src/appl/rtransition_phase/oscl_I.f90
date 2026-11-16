      module oscl_i
      implicit none (type, external)
      interface
      subroutine oscl (name)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      end subroutine oscl
      end interface
      end module oscl_i
