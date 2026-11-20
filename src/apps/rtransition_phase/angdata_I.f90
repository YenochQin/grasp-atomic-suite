      module angdata_i
      implicit none (type, external)
      interface
      subroutine angdata (name, avail, jkp, transition_angular_unit)
      implicit none (type, external)
      character (len = 24), intent(inout) :: name(2)
      logical, intent(out) :: avail
      integer :: jkp
      integer, intent(out) :: transition_angular_unit
      end subroutine angdata
      end interface
      end module angdata_i
