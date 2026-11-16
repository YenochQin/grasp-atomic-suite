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
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only:  max_radial_points
      use memory_man
      use debug_c
      use def_c, only: speed_of_light_au, nuclear_charge
      use grid_c
      use npar_c
      use orb_c
      use wave_c, only: orbital_origin_parameters, radial_large_components, radial_small_components
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpq_i
      use orthsc_i
      use mpiu_i, only: stopmpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: ierror
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npy, naky, my
!     integer :: ierr
      real(real64) :: con, fkk, ey, pzy, dnorm
      real(real64), pointer :: pa(:), qa(:), ra(:)
!-----------------------------------------------
!
!-----------------------------------------------------------------------
!
!   Allocate storage to orbital arrays
!
      call alloc(radial_large_components, max_radial_points, orbital_count, 'PF', 'LODRWFMPI')
      call alloc(radial_small_components, max_radial_points, orbital_count, 'QF', 'LODRWFMPI')

!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      con = nuclear_charge / speed_of_light_au
      con = con * con

      do j = 1, orbital_count

         do i = 1, radial_point_count
            radial_large_components(i,j) = 0._real64
            radial_small_components(i,j) = 0._real64
         enddo

         orbital_energy_parameters(j) = -1._real64
         k = abs (orbital_kappa(j))
         if (nuclear_model_parameter_count .gt. 0) then
            orbital_origin_exponents(j) = dble (k)
         elseif (nuclear_model_parameter_count .eq. 0) then
            fkk = dble (k*k)
            if (fkk .ge. con) then
               orbital_origin_exponents(j) = sqrt (fkk-con)
            else
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               call stopmpi ('lodrwfmpi: Inappropriate gamma', mpi_rank)
            endif
         endif
      enddo
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      if (ldbpr(3) .and. mpi_rank .eq. 0) write (99,300)
      nwin = 0
    3 continue

      if (mpi_rank .eq. 0) read (23,iostat = ios) npy,naky,ey,my
      call mpi_bcast (ios,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (npy,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (naky,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (my,1,mpi_integer,0,mpi_comm_world,mpi_error_code)
      call mpi_bcast (ey,1,mpi_double_precision,0,mpi_comm_world,mpi_error_code)

      if (ios .eq. 0) then
         call alloc (pa,my, 'PA', 'LODRWFMPI')
         call alloc (qa,my, 'QA', 'LODRWFMPI')
         call alloc (ra,my, 'RA', 'LODRWFMPI')

         if (mpi_rank .eq. 0) then
            read (23) pzy,(pa(i),i = 1,my),(qa(i),i =1 ,my)
            read (23) (ra(i),i = 1,my)
         endif

         call mpi_bcast (pzy,1,mpi_double_precision,0, &
                               mpi_comm_world,mpi_error_code)
         call mpi_bcast (pa,my,mpi_double_precision,0, &
                               mpi_comm_world,mpi_error_code)
         call mpi_bcast (qa,my,mpi_double_precision,0, &
                               mpi_comm_world,mpi_error_code)
         call mpi_bcast (ra,my,mpi_double_precision,0, &
                               mpi_comm_world,mpi_error_code)

         do j = 1, orbital_count
            if (.not.(orbital_energy_parameters(j)<0.0_real64 .and. npy==orbital_principal_n(j) .and. naky==orbital_kappa(j))) cycle
               orbital_origin_parameters(j) = pzy
               orbital_energy_parameters(j) = ey
               call intrpq (pa, qa, my, ra, j, dnorm)
               if (ldbpr(3) .and. mpi_rank .eq. 0) &
                     write (99,301) orbital_principal_n(j), orbital_symmetry_labels(j), orbital_energy_parameters(j), dnorm
               nwin = nwin + 1
         enddo

         call dalloc (pa, 'PA', 'LODRWFMPI')
         call dalloc (qa, 'QA', 'LODRWFMPI')
         call dalloc (ra, 'RA', 'LODRWFMPI')

         goto 3
      endif
      if (ldbpr(3) .and. mpi_rank .eq. 0) &
            write (99,*) ' orbitals renormalised;'
!
!   Return with an error code if all orbitals are not known
!
      if (nwin .lt. orbital_count) then
         ierror = 1
         return
      endif
!
!   Schmidt orthogonalise the orbitals
!
      call orthsc

      if (ldbpr(3) .and. mpi_rank .eq. 0) &
            write (99,*) 'orbitals orthogonalised and renormalised;'
      ierror = 0

  300 format (/'From SUBROUTINE LODRWFmpi:' &
              /' Orbital',8x,'Eigenvalue',19x,'Norm')
  301 format (2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)

      return
      end subroutine lodrwfmpi
