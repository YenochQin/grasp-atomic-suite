!***********************************************************************
!                                                                      *
      subroutine lodrwfi(name, ntestg)
      use calculation_state_files, only: orbital_input_path
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Per Jonsson                              June 1996      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_orbitals, max_radial_points
      use array_memory_management, only: alloc, dalloc
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use atomic_calculation_state,           only: nuclear_charge, speed_of_light_au
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use biorthogonal_subshell_layout,         only: kamax, nshlii, nshlpii
      use radial_wavefunctions,          only: pfii, qfii, pzii, mfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpqi_i, only: intrpqi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ntestg
      character , intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  :: nak(max_orbitals), np(max_orbitals)
      integer :: ntestl, ntest, j, k, i, nwin, ios, npy, naky, my, jj, kk, jjj&
         , kkk, ierr
      real(real64) :: e(max_orbitals), gama(max_orbitals)
      real(real64) :: con, fkk, ey, pzy, dnorm
      real(real64), pointer :: pa(:), qa(:), ra(:)
      character  :: nh(max_orbitals)*2
      character :: g92rwf*6
!-----------------------------------------------
!
!  Common relevant for the initial state
!
!
      ntestl = 0
      ntest = max0(ntestl,ntestg)
      ntest = 0

!
!   Write entry message
!
      write (6, *) "Loading Radial WaveFunction File for initial state..."
!
!   Open the radial wave function file
!
      open(unit=21, file=orbital_input_path(name), form="UNFORMATTED", status="OLD", &
         position="asis")
!
!   Save NAK, NP and NH
!
      nak(:nwii) = nakii(:nwii)
      np(:nwii) = npii(:nwii)
      nh(:nwii) = nhii(:nwii)
!
!   Allocate storage to orbital arrays
!
      call alloc (pfii, max_radial_points,nwii, "PFII", "LODRWII")
      call alloc (qfii, max_radial_points,nwii, "QFII", "LODRWII")

!
      con = nuclear_charge/speed_of_light_au
      con = con*con
!
      do j = 1, nwii
         pfii(:max_radial_points,j) = 0.0_real64
         qfii(:max_radial_points,j) = 0.0_real64
!
         k = abs(nak(j))
         if (nuclear_model_parameter_count /= 0) cycle
         fkk = dble(k*k)
         if (fkk >= con) then
            gama(j) = sqrt(fkk - con)
         else
            write (6, *) "LODRWF: Imaginary gamma parameter"
            write (6, *) " for ", np(j), nh(j), " orbital; the"
            write (6, *) " point model for the nucleus"
            write (6, *) " is inappropriate for Z > ", speed_of_light_au, "."
            stop
         end if
!
      end do
!
!   Read orbital information from Read Orbitals File;
!
      nwin = 0
      read (21, iostat=ios) g92rwf
      if (ios/=0 .or. g92rwf/="G92RWF") then
         write (6, *) "This is not a Radial WaveFunction File;"
         close(21)
      end if

      if (ntest >= 100) then
         write (*, *) "******************"
         write (*, *) " Entering lodrwfi"
         write (*, *) "******************"
      end if

    3 continue
      read (21, iostat=ios) npy, naky, ey, my
      if (ios == 0) then
         call alloc (pa, my, "PA", "LODRWFII")
         call alloc (qa, my, "QA", "LODRWFII")
         call alloc (ra, my, "RA", "LODRWFII")
         read (21) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
         read (21) (ra(i),i=1,my)
!
!    Orbital order as defined in kapdata
!
         jj = 0
         do k = 1, kamax
            if (k > 1) jj = nshlii(k-1) + jj
            do j = 1, nshlii(k)
               kk = nshlpii(k,j)
               if (npy/=np(kk) .or. naky/=nak(kk)) cycle
               jjj = jj + j
               pzii(jjj) = pzy
               eii(jjj) = ey
               nakii(jjj) = nak(kk)
               npii(jjj) = np(kk)
               nhii(jjj) = nh(kk)
               gamaii(jjj) = gama(kk)
               call intrpqi (pa, qa, my, ra, jjj, dnorm)
               if (ntest >= 100) write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') npii(jjj), nhii(jjj), eii(jjj)&
                  , dnorm
               if (ntest > 1000) then
                  write (*, *) "PF              QF             RA"
                  do kkk = 1, mfii(jjj)
                     write (*, *) pfii(kkk,jjj), qfii(kkk,jjj), ra(kkk)
                  end do
               end if
               nwin = nwin + 1
            end do
         end do
         call dalloc (pa, "PA", "LODRWII")
         call dalloc (qa, "QA", "LODRWII")
         call dalloc (ra, "RA", "LODRWII")
         goto 3
      end if
      if (ldbpr(3)) write (99, *) " orbitals renormalised;"
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nwii) then
         write (6, *) "LODRWF: All required orbitals not"
         write (6, *) " found."
         ierr = 1
         goto 5
      end if
!
      write (6, *) " ... load complete;"
!
    5 continue
      close(21)
      if (ntest >= 100) then
         write (*, *) "Sorted order should be the same as from kapdat"
         do j = 1, nwii
            write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') npii(j), nhii(j), eii(j), dnorm
         end do
         write (*, *)
         write (*, *) "*****************"
         write (*, *) " Leaving lodrwfi"
         write (*, *) "*****************"
      end if

      return
!
      return
!
      end subroutine lodrwfi
