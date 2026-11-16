      module iniestsd_i
      implicit none (type, external)
      interface
      subroutine iniestsd (nmax,ncf, myid,nprocs,niv,basis,imcdf,eav)
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      integer, intent(in) :: nmax, ncf, myid, nprocs, niv, imcdf
      real(real64), intent(in) :: eav
      real(real64), intent(inout):: basis(*)
      end subroutine
      end interface
      end module
