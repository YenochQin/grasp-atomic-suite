      module mctout_i
      implicit none (type, external)
      interface
      subroutine mctout (operator_parity, transition_type_index, name)
      implicit none (type, external)
      integer, intent(out) :: operator_parity
      integer, intent(in) :: transition_type_index
      character (len = 24) :: name(2)
      end subroutine mctout
      end interface
      end module mctout_i
