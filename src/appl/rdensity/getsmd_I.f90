      module getsmd_i
      implicit none (type, external)
      interface
      subroutine getsmd (name)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      end subroutine getsmd
      end interface
      end module getsmd_i
