      module setsum_i
      implicit none (type, external)
      interface
      subroutine setsum(defnam)
      implicit none (type, external)
      character(len=*), intent(in) :: defnam
      end subroutine setsum
      end interface
      end module setsum_i
