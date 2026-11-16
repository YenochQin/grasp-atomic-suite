      module sbdat1_c
      use vast_kind_param, only: double
      use parameter_def,   only: max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  06:25:32  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: nlmax = 20
      integer, dimension(nlmax,nlmax) :: nshlp
      integer, dimension(nlmax,max_orbitals) :: nshlpp
      end module sbdat1_c
