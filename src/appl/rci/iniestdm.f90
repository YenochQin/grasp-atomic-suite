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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
      use eigv_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dspevx_i
      use suite_linear_algebra, only: dcopy, dspevx
      use dinit_i, only: dinit
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in):: nmax, ncf, niv
      real(double), dimension(*) :: basis, hmx
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------------------------------
      integer :: myid, nprocs, info, ierr, ns, i, nfound, j
      real(double), dimension(:), pointer :: ap, eigval, vec, work
      integer, dimension(:), pointer :: iwork, ifail
!-----------------------------------------------------------------------

      myid = 0
      nprocs = 1

      ns = min (nmax, ncf)

      call alloc (ap, ns*(ns+1)/2, 'AP', 'INIESTDM')
      call dinit (ns*(ns+1)/2, 0.d0, ap, 1)

!  Get the upper left sub-matrix

      do i = 1, (ns*(ns+1))/2
         ap(i) = hmx(i)
      enddo

      call alloc (eigval,   ns, 'EIGVAL', 'INIESTDM')
      call alloc (vec,   ns*niv, 'VEC', 'INIESTDM')
      call alloc (work,  8*ns,  'WORK', 'INIESTDM')
      call alloc (iwork, 5*ns,  'IWORK', 'INIESTDM')
      call alloc (ifail,     ns,  'IFAIL', 'INIESTDM')

      call dspevx ('Vectors also','In a range','Upper triangular',      &
     &          ns,ap,-1.0d0,-1.0d0,1,niv,0.d0,                               &
     &          nfound,eigval,vec,ns,work,iwork,ifail,info)
           ierr = -abs (info)

!  Build the Basis.

      call dinit (ncf*niv, 0.d0, basis, 1)

!  scatter the vectors

      do j = 1, niv
            call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)),1, basis(ncf*(j-1)+1), 1)
      enddo

      call dcopy (niv, eigval, 1, basis(niv*ncf+1), 1)

      call dalloc (ap, 'AP', 'INIESTDM')
      call dalloc (eigval, 'EIGVAL', 'INIESTDM')
      call dalloc (vec, 'VEC', 'INIESTDM')
      call dalloc (work, 'WORK', 'INIESTDM')
      call dalloc (iwork, 'IWORK', 'INIESTDM')
      call dalloc (ifail, 'IFAIL', 'INIESTDM')

      return
      end
