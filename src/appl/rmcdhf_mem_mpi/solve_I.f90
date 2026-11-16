      module solve_i
      implicit none (type, external)
      interface
      subroutine solve (j, fail, inv, jp, nnp)
      implicit none (type, external)
      integer, intent(in) :: j
      logical, intent(out) :: fail
      integer, intent(inout) :: inv
      integer, intent(out) :: jp
      integer, intent(out) :: nnp
      end subroutine solve
      end interface
      end module solve_i
