      module hydrogenic_screening
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_orbitals
      implicit none (type, external)
      private :: real64
      real(real64) :: sigma(max_orbitals)
      end module hydrogenic_screening
