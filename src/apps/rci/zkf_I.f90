      module zkf_i
      implicit none (type, external)
      interface
      subroutine zkf (k, i, j)
      implicit none (type, external)
      integer, intent(in) :: k
      integer :: i
!VAST...Dummy argument I is not referenced in this routine.
      integer :: j
      end subroutine zkf
      end interface
      end module zkf_i
