      module scf_iteration_state
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_orbitals
      implicit none (type, external)
      private :: real64
      real(real64) :: ucf(max_orbitals)
      integer :: method(max_orbitals)
      real(real64) :: scnsty(max_orbitals)
      real(real64) :: epsmin, epsmax, emin, emax, zinf
      integer :: ndcof, nxcof, nycof, nddim, nxdim, nydim
      real(real64), pointer :: da(:), xa(:), ya(:)
      integer, pointer :: nda(:), nxa(:), nya(:)
      end module scf_iteration_state
