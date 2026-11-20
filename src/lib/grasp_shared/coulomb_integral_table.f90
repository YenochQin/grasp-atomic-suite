      module coulomb_integral_table
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: ndcoea, ncoei
      integer, pointer :: indoei(:)
      real(real64), pointer :: valoei(:)
      logical :: frstco
      end module coulomb_integral_table
