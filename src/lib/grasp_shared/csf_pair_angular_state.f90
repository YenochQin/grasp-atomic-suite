! Holds occupations, coupling data, and the peel/core orbital lists for the current CSF pair.
      module csf_pair_angular_state
      use grasp_parameters,   only:  max_orbitals
      implicit none (type, external)
      integer :: nq1(max_orbitals), nq2(max_orbitals)
      integer :: jjc1(max_orbitals), jjc2(max_orbitals)
      integer :: jjq1(3, max_orbitals), jjq2(3, max_orbitals)
      integer :: jlist(max_orbitals), klist(max_orbitals)
      integer :: npeel, ncore
      end module csf_pair_angular_state
