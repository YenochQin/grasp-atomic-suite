      module fixed_orbital_flags
      use grasp_parameters,    only: max_orbitals
      implicit none (type, external)
      integer :: fixed_orbital_count
      logical :: is_orbital_fixed(max_orbitals)
      end module fixed_orbital_flags
