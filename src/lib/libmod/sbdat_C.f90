      module sbdat_c
      use parameter_def,   only: max_orbitals
      implicit none (type, external)
!...Created by Pacific-Sierra Research 77to90  4.3E  06:27:59  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer, parameter :: nlmax = 20
      integer :: nakinvii(max_orbitals)
      integer :: nshlii(nlmax)
      integer :: nshlpii(nlmax, nlmax)
      integer :: nakinvff(max_orbitals)
      integer :: nshlff(nlmax)
      integer :: nshlpff(nlmax, nlmax)
      integer :: nshlppii(nlmax, max_orbitals), nshlppff(nlmax, max_orbitals)
      integer :: ninii(nlmax), ninff(nlmax), ikappa(nlmax)
      integer :: kamax
      end module sbdat_c
