      module iniestdm_i
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine iniestdm (nmax, ncf, niv, basis, hmx)
      use vast_kind_param,only: double
      integer, intent(in):: nmax, ncf, niv
      real(double), dimension(*) :: basis, hmx
      end subroutine  iniestdm
      end interface
      end module
