      module orbital_lagrange_multipliers
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: nec
      real(real64), pointer :: ecv(:)
      integer, pointer :: iecc(:)
      end module orbital_lagrange_multipliers
