      module fixd_c
      use vast_kind_param,  only: double
      use parameter_def,    only: max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  06:37:37  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: fixed_orbital_count
      logical, dimension(max_orbitals) :: is_orbital_fixed
      end module fixd_c
