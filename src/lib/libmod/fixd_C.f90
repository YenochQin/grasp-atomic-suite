      module fixd_c
      use, intrinsic :: iso_fortran_env, only: double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:37:37  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: fixed_orbital_count
      logical :: is_orbital_fixed(max_orbitals)
      end module fixd_c
