! Holds energy shifts, eigenvalues, and CSF mixing coefficients for the active, final, and initial states.
      module eigenstate_expansions
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      real(real64) :: hamiltonian_energy_shift, eavff, eavii
      real(real64), pointer :: shifted_eigenvalues(:), evalff(:), evalii(:)
      real(real64), pointer :: csf_mixing_coefficients(:), evecff(:), evecii(:)
      end module eigenstate_expansions
