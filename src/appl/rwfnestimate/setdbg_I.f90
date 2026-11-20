      module setdbg_i
      implicit none (type, external)
      interface
      subroutine setdbg(defnam)
      implicit none (type, external)
      character(len=*), intent(in) :: defnam
      end subroutine setdbg
      end interface
      end module setdbg_i
