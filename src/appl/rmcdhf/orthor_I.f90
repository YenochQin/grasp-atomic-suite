      module orthor_i
      implicit none (type, external)
      interface
      subroutine orthor (j, inv)
      implicit none (type, external)
      integer, intent(in) :: j
      integer, intent(inout) :: inv
      end subroutine orthor
      end interface
      end module orthor_i
