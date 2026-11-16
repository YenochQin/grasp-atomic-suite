      module ris_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      implicit none (type, external)
      private :: real64
!...Translated by Gediminas Gaigalas 11/18/19
      real(real64), pointer :: dens1(:), dens2(:), dens3(:),  &
                                             dens4(:), dens5(:), dens6(:),  &
                                             dens7(:), smsc1(:), smsc2(:)
      end module ris_c
