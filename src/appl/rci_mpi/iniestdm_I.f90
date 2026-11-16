      module iniestdm_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine iniestdm (nmax, ncf, niv, basis, hmx)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in):: nmax, ncf, niv
      real(real64) :: basis(*), hmx(*)
      end subroutine  iniestdm
      end interface
      end module
