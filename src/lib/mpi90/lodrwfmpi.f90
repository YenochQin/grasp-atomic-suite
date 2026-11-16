!***********************************************************************
      subroutine lodrwfmpi (ierror)
!     IMPLICIT REAL*8          (A-H, O-Z)

!   This subroutine loads radial wavefunctions from the .rwf  file
!   and performs some related setup. It does not handle any error
!   so the caller has to do it.
!       ierror=0, normal
!       ierror=1, error
!
!   Used by rcimpivu, rscfmpi, rcimpi
!
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
!   MPI version by Xinghong He            Last revision: 27 May 1997   *
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only:  double
      use parameter_def,   only:  nnnp
      use memory_man
      use debug_c
      use def_c, only: c, z
      use grid_c
      use npar_c
      use orb_c
      use wave_c, only: pz, pf, qf
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpq_i
      use orthsc_i
      use mpiu_i, only: stopmpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: ierror
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npy, naky, my
!     integer :: ierr
      real(double) :: con, fkk, ey, pzy, dnorm
      real(double), dimension(:), pointer :: pa, qa, ra
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!
!   Allocate storage to orbital arrays
!
      call alloc(pf, nnnp, nw, 'PF', 'LODRWFMPI')
      call alloc(qf, nnnp, nw, 'QF', 'LODRWFMPI')

!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      con = z / c
      con = con * con

      do j = 1, nw

         do i = 1, n
            pf(i,j) = 0.d0
            qf(i,j) = 0.d0
         enddo

         e(j) = -1.d0
         k = abs (nak(j))
         if (nparm .gt. 0) then
            gama(j) = dble (k)
         elseif (nparm .eq. 0) then
            fkk = dble (k*k)
            if (fkk .ge. con) then
               gama(j) = sqrt (fkk-con)
            else
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               call stopmpi ('lodrwfmpi: Inappropriate gamma', myid)
            endif
         endif
      enddo
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      if (ldbpr(3) .and. myid .eq. 0) write (99,300)
      nwin = 0
    3 continue

      if (myid .eq. 0) read (23,iostat = ios) npy,naky,ey,my
      call mpi_bcast (ios,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (npy,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (naky,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (my,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast (ey,1,mpi_double_precision,0,mpi_comm_world,ierr)

      if (ios .eq. 0) then
         call alloc (pa,my, 'PA', 'LODRWFMPI')
         call alloc (qa,my, 'QA', 'LODRWFMPI')
         call alloc (ra,my, 'RA', 'LODRWFMPI')

         if (myid .eq. 0) then
            read (23) pzy,(pa(i),i = 1,my),(qa(i),i =1 ,my)
            read (23) (ra(i),i = 1,my)
         endif

         call mpi_bcast (pzy,1,mpi_double_precision,0, &
                               mpi_comm_world,ierr)
         call mpi_bcast (pa,my,mpi_double_precision,0, &
                               mpi_comm_world,ierr)
         call mpi_bcast (qa,my,mpi_double_precision,0, &
                               mpi_comm_world,ierr)
         call mpi_bcast (ra,my,mpi_double_precision,0, &
                               mpi_comm_world,ierr)

         do j = 1, nw
            if (.not.(e(j)<0.0d00 .and. npy==np(j) .and. naky==nak(j))) cycle
               pz(j) = pzy
               e(j) = ey
               call intrpq (pa, qa, my, ra, j, dnorm)
               if (ldbpr(3) .and. myid .eq. 0) &
                     write (99,301) np(j), nh(j), e(j), dnorm
               nwin = nwin + 1
         enddo

         call dalloc (pa, 'PA', 'LODRWFMPI')
         call dalloc (qa, 'QA', 'LODRWFMPI')
         call dalloc (ra, 'RA', 'LODRWFMPI')

         goto 3
      endif
      if (ldbpr(3) .and. myid .eq. 0) &
            write (99,*) ' orbitals renormalised;'
!
!   Return with an error code if all orbitals are not known
!
      if (nwin .lt. nw) then
         ierror = 1
         return
      endif
!
!   Schmidt orthogonalise the orbitals
!
      call orthsc

      if (ldbpr(3) .and. myid .eq. 0) &
            write (99,*) 'orbitals orthogonalised and renormalised;'
      ierror = 0

  300 format (/'From SUBROUTINE LODRWFmpi:' &
              /' Orbital',8x,'Eigenvalue',19x,'Norm')
  301 format (2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)

      return
      end subroutine lodrwfmpi
