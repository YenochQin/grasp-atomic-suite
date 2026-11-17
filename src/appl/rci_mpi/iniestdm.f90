!************************************************************************
!
      subroutine iniestdm (nmax, ncf, niv, basis, hmx)

!     Routine for providing initial estimates from the diagonal
!     of the matrix. This way was used by Dvdson in atomic structure
!     calculations. It should be used to obtain estimates when nothing
!     else is available.
!
!   Block version by Xinghong He            Last revision: 18 Jun 1998
!
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man, only: alloc, dalloc, double
      use eigv_c, only: double, hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use mpi_c, only: double, mpi_error_code, mpi_process_count, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dspevx_i, only: dspevx
      use suite_linear_algebra, only: dcopy, dspevx
      use dinit_i, only: dinit
      use suite_linear_algebra, only: slamch
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in):: nmax, ncf, niv
      real(real64) :: basis(*), hmx(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      integer :: info, ns, i, nfound, j, ihmx, iap
      real(real64), pointer :: ap(:), eigval(:), vec(:), work(:)
      integer, pointer :: iwork(:), ifail(:)
      real(real64) :: abstol
!-----------------------------------------------------------------------
      ns = min (nmax, ncf)

      call alloc (ap, ns*(ns+1)/2, "AP", "INIESTDM")
      call dinit (ns*(ns+1)/2, 0.0_real64, ap, 1)

!  Get the upper left sub-matrix
!  ap is global, hmx is local !!!

      ihmx = 0
      do i = mpi_rank + 1, ns, mpi_process_count
         iap = i*(i-1)/2
         do j = 1, i
            ap(j+iap) = hmx(j+ihmx)
         end do
         ihmx = ihmx + i
      end do

      call alloc (eigval,   ns, "EIGVAL", "INIESTDM")
      call alloc (vec,   ns*niv, "VEC", "INIESTDM")
      call alloc (work,  8*ns,  "WORK", "INIESTDM")
      call alloc (iwork, 5*ns,  "IWORK", "INIESTDM")
      call alloc (ifail,     ns,  "IFAIL", "INIESTDM")

!      CALL DSPEVX ('Vectors also','In a range','Upper triangular',     &
!     &          NS,AP,-1.0D0,-1.0D0,1,NIV,0.d0,                        &
!     &          NFOUND,EIGVAL,VEC,NS,work,iwork,IFAIL,INFO)
      abstol = 2*slamch("S")
      call dspevx ("V","I","U",                                        &
                ns,ap,-1.0_real64,-1.0_real64,1,niv,abstol,                            &
                nfound,eigval,vec,ns,work,iwork,ifail,info)
      mpi_error_code = -abs (info)

!  Build the Basis.

      call dinit (ncf*niv, 0.0_real64, basis, 1)

!  scatter the vectors

      do j = 1, niv
            call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)),1, basis(ncf*(j-1)+1), 1)
      end do

      call dcopy (niv, eigval, 1, basis(niv*ncf+1), 1)

      call dalloc (ap, "AP", "INIESTDM")
      call dalloc (eigval, "EIGVAL", "INIESTDM")
      call dalloc (vec, "VEC", "INIESTDM")
      call dalloc (work, "WORK", "INIESTDM")
      call dalloc (iwork, "IWORK", "INIESTDM")
      call dalloc (ifail, "IFAIL", "INIESTDM")

      return
      end subroutine iniestdm
