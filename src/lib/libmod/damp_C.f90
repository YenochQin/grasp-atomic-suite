      module damp_c
!     Arrays associated with damping solutions in the SCF process
!     ... ODAMP:  array for damping radial functions
!     ... CDAMP:  array for damping expansion coefficients
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
      private :: real64
      real(real64) :: orbital_damping_factors(max_orbitals)
      real(real64), pointer :: mixing_damping_factors(:)
      end module damp_c
