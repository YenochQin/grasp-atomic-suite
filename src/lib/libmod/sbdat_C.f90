      module sbdat_c
      use parameter_def,   only: max_orbitals
!...Created by Pacific-Sierra Research 77to90  4.3E  06:27:59  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: nlmax = 20
      integer, dimension(max_orbitals) :: nakinvii
      integer, dimension(nlmax) :: nshlii
      integer, dimension(nlmax,nlmax) :: nshlpii
      integer, dimension(max_orbitals) :: nakinvff
      integer, dimension(nlmax) :: nshlff
      integer, dimension(nlmax,nlmax) :: nshlpff
      integer, dimension(nlmax,max_orbitals) :: nshlppii, nshlppff
      integer, dimension(nlmax) :: ninii, ninff, ikappa
      integer :: kamax
      end module sbdat_c
