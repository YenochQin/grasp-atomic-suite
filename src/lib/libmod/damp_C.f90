      module damp_c
!     Arrays associated with damping solutions in the SCF process
!     ... ODAMP:  array for damping radial functions
!     ... CDAMP:  array for damping expansion coefficients
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      use vast_kind_param,  only: double
      use parameter_def,    only: nnnw
      real(double), dimension(nnnw) :: odamp
      real(double), dimension(:), pointer :: cdamp
      end module damp_c
