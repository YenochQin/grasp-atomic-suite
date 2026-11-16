      module mctin_i
      implicit none (type, external)
      interface
      subroutine mctin (iopar, jkp, name)
      implicit none (type, external)
      integer :: iopar
!VAST...Dummy argument IOPAR is not referenced in this routine.
      integer, intent(in) :: jkp
      character (len = 24) :: name(2)
      end subroutine mctin
      end interface
      end module mctin_i
