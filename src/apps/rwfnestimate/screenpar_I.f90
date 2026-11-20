      module screenpar_i
      implicit none (type, external)
      interface
      subroutine screenpar (core_subshell_count)
      use grasp_parameters, only: max_orbitals
      implicit none (type, external)
      integer, intent(in) :: core_subshell_count
      end subroutine screenpar
      end interface
      end module screenpar_i
