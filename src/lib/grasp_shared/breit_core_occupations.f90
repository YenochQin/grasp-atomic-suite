! Stores core occupations used by the Breit matrix-element construction.
      module breit_core_occupations
      use grasp_parameters,   only:  max_orbitals
      implicit none (type, external)
      integer :: icore(max_orbitals)
      end module breit_core_occupations
