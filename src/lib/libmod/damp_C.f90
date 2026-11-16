      module damp_c
!     Arrays associated with damping solutions in the SCF process
!     ... ODAMP:  array for damping radial functions
!     ... CDAMP:  array for damping expansion coefficients
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use vast_kind_param,  only: double
      use parameter_def,    only: max_orbitals
      real(double), dimension(max_orbitals) :: orbital_damping_factors
      real(double), dimension(:), pointer :: mixing_damping_factors
      end module damp_c
