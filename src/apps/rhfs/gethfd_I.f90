      module gethfd_i
      implicit none (type, external)
      interface
      subroutine gethfd (name)
      implicit none (type, external)
      character (len = 24), intent(in) :: name
      end subroutine gethfd
      end interface
      end module gethfd_i
