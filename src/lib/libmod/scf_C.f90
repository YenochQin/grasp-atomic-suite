      module scf_c
      use, intrinsic :: iso_fortran_env, only: real64, double => real64
      use parameter_def,    only: max_orbitals
      implicit none (type, external)
      private :: real64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      real(real64) :: ucf(max_orbitals)
      integer :: method(max_orbitals)
      real(real64) :: scnsty(max_orbitals)
      real(real64) :: epsmin, epsmax, emin, emax, zinf
      integer :: ndcof, nxcof, nycof, nddim, nxdim, nydim
      real(real64), pointer :: da(:), xa(:), ya(:)
      integer, pointer :: nda(:), nxa(:), nya(:)
      end module scf_c
