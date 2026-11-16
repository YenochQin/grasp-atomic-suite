!***********************************************************************
!
      subroutine iniestsd (nmax, ncf, myid, nprocs,                    &
          niv, basis, imcdf, eav)
!
!    Routine for providing initial estimates from upper left corner
!       of the matrix. It is exact (not estimates) if ncf <= nmax which
!    is currently set to be 400 in the calling routine.
!
!    Matrix is sparse and on the disk
!
!   Block version by Xinghong He            Last revision: 14 Dec 1998
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use memory_man
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
       integer, intent(in) :: nmax, ncf, myid, nprocs, niv, imcdf
       real(double), intent(in) :: eav
       real(double), dimension(*), intent(inout):: basis
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
       real(double), dimension(:), pointer :: ap, eigval, vec, work, hmx
       integer, dimension(:), pointer :: iwork, ifail, irow
       integer :: in, ns, ncfdum, iccutdum, myiddum, nprocsdum, info, &
            ierr, j, joff, nelc, ir, nfound
       real(double) :: elsto
!-----------------------------------------------------------------------
      ns = min (nmax, ncf)

      call alloc (ap, (ns*(ns+1))/2, 'AP', 'INIESTSD')
      call dinit ((ns*(ns+1))/2, 0.d0, ap, 1)

!**** separate upper left block of size NS*NS

      call alloc (hmx, ncf, 'HMX', 'INIESTSD')
      call alloc (irow, ncf, 'IROW', 'INIESTSD')
      read (imcdf) ncfdum, iccutdum, myiddum, nprocsdum
      if (myid .eq. 0) print *, 'iniestsd ...........'
      if (ncf .ne. ncfdum .or.  myid .ne. myiddum                       &
                    .or. nprocsdum .ne. nprocs)                         &
         stop 'iniestsd: ncf read wrong'

!          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!          if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n.

      do j = myid + 1, ns, nprocs
         joff = (j*(j-1))/2
         read (imcdf) nelc,elsto,(hmx(ir),ir=1,nelc),                   &
                                (irow(ir),ir=1,nelc)
         hmx(nelc) = hmx(nelc) - eav ! Shift the diagonal
         do ir = 1, nelc
            ap(irow(ir) + joff) = hmx(ir)
         enddo
      enddo

! Let each node have a complete copy of ap

!      CALL gdsummpi (ap, (NS*(NS+1))/2)

! To be in step with other cases, go through the whole block.
!
! This is not necessary since currently the file pointer is moved
! to the absolute position of the .res files which is always counted
! from the begining of the .res files of each node. Besides, the
! following segment seems not working properly for the last block.
! Xinghong He 98-12-14

      !mylast = j - nprocs
      !DO j = mylast, ncf, nprocs
      !   READ (imcdf)
      !ENDDO

      call dalloc (hmx,'HMX', 'INIESTSD')
      call dalloc (irow, 'HMX', 'INIESTSD')

      call alloc (eigval,ns,'EIGVAL','INIESTSD')
      call alloc (vec,ns*niv,'VEC','INIESTSD')
      call alloc (work,8*ns,'WORK','INIESTSD')
      call alloc (iwork,5*ns,'IWORK','INIESTSD' )
      call alloc (ifail,ns, 'IFAIL','INIESTSD' )

      call dspevx ('Vectors also','In a range','Upper triangular',     &
                ns,ap,-1.0d0,-1.0d0,1,niv,0.d0,                        &
                nfound,eigval,vec,ns,work,iwork,ifail,info)
      ierr = -abs (info)

!******************************************************************

!       ..Build the Basis.

      call dinit (ncf*niv, 0.d0, basis, 1)
!       ...scatter the vectors
      do j = 1, niv
         call dcopy (ns, vec(ns*(j-1)+1:ns*(j-1)+(ns)),1, basis(ncf*(j-1)+1), 1)
      enddo
      call dcopy (niv, eigval,1,basis(niv*ncf+1),1)

      call dalloc (ap, 'AP', 'INIESTSD')
      call dalloc (eigval,'EIGVAL', 'INIESTSD')
      call dalloc (vec, 'VEC', 'INIESTSD')
      call dalloc (work, 'WORK', 'INIESTSD')
      call dalloc (iwork, 'IWORK', 'INIESTSD')
      call dalloc (ifail, 'IFAIL', 'INIESTSD')

      return
      end subroutine iniestsd
