      module iniestsd_i
      implicit none (type, external)
      interface
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      subroutine iniestsd (nmax,ncf, myid,nprocs,niv,basis,imcdf,eav)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nmax, ncf, myid, nprocs, niv, imcdf
      real(real64), intent(in) :: eav
      real(real64), intent(inout):: basis(*)
      end subroutine
      end interface
      end module
