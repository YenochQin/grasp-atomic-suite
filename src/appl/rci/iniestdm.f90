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
      use memory_man
      use eigv_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dspevx_i
      use suite_linear_algebra, only: dcopy, dspevx
      use dinit_i, only: dinit
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in):: nmax, ncf, niv
      real(real64) :: basis(*), hmx(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      integer :: myid, nprocs, info, ierr, ns, i, nfound, j
      real(real64), pointer :: ap(:), eigval(:), vec(:), work(:)
      integer, pointer :: iwork(:), ifail(:)
!-----------------------------------------------------------------------

      myid = 0
      nprocs = 1

      ns = min (nmax, ncf)

      call alloc (ap, ns*(ns+1)/2, "AP", "INIESTDM")
      call dinit (ns*(ns+1)/2, 0.0_real64, ap, 1)

!  Get the upper left sub-matrix

      do i = 1, (ns*(ns+1))/2
         ap(i) = hmx(i)
      end do

      call alloc (eigval,   ns, "EIGVAL", "INIESTDM")
      call alloc (vec,   ns*niv, "VEC", "INIESTDM")
      call alloc (work,  8*ns,  "WORK", "INIESTDM")
      call alloc (iwork, 5*ns,  "IWORK", "INIESTDM")
      call alloc (ifail,     ns,  "IFAIL", "INIESTDM")

      call dspevx ("Vectors also","In a range","Upper triangular",      &
     &          ns,ap,-1.0_real64,-1.0_real64,1,niv,0.0_real64,                               &
     &          nfound,eigval,vec,ns,work,iwork,ifail,info)
           ierr = -abs (info)

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
