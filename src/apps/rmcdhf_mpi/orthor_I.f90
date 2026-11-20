      module orthor_i
      implicit none (type, external)
      interface
      subroutine orthor (orbital_index, orbital_sign_reversal_count)
      implicit none (type, external)
      integer, intent(in) :: orbital_index
      integer, intent(inout) :: orbital_sign_reversal_count
      end subroutine orthor
      end interface
      end module orthor_i
