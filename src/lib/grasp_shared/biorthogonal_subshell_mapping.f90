      module biorthogonal_subshell_mapping
      use grasp_parameters,   only: max_orbitals
      implicit none (type, external)
      integer, parameter :: nlmax = 20
      integer :: nshlp(nlmax, nlmax)
      integer :: nshlpp(nlmax, max_orbitals)
      end module biorthogonal_subshell_mapping
