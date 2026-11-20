      module kinetic_energy_integral_table
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: ndkea, nkei
      integer, pointer :: indkei(:)
      real(real64), pointer :: valkei(:)
      logical :: frstki
      end module kinetic_energy_integral_table
