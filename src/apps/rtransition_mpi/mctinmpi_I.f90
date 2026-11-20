      module mctin_i
      implicit none (type, external)
      interface
      subroutine mctin (operator_parity, transition_type_index, name, transition_angular_unit)
      implicit none (type, external)
      integer, intent(in) :: transition_angular_unit
      integer :: operator_parity
!VAST...Dummy argument IOPAR is not referenced in this routine.
      integer, intent(in) :: transition_type_index
      character (len = 24) :: name(2)
      end subroutine mctin
      end interface
      end module mctin_i
