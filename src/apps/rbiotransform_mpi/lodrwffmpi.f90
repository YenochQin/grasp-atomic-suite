!***********************************************************************
!                                                                      *
      subroutine lodrwff(state_name, ntestg)
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
      use biorthogonal_subshell_layout,         only: kamax, nshlff, nshlpff
      use radial_wavefunctions,          only: pfff, qfff, pzff, mfff
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpqf_i, only: intrpqf
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ntestg
      character, intent(in) :: state_name*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  :: nak(max_orbitals), np(max_orbitals)
      integer :: ntestl, ntest, j, k, i, nwin, io_status, npy, naky, my, jj, kk, jjj&
         , kkk, ierr
      real(real64) :: e(max_orbitals), gama(max_orbitals)
      real(real64) :: con, fkk, ey, dnorm, pzy
      real(real64), pointer :: pa(:), qa(:), ra(:)
      character  :: nh(max_orbitals)*2
      character :: radial_file_signature*6
!-----------------------------------------------
!
!   Common relevant for the final state
!
!
      ntestl = 0
      ntest = max0(ntestl,ntestg)
      ntest = 0

!
!   Write entry message
!
      write (6, *) "Loading Radial WaveFunction File for final state..."
!
!   Open the radial wave function file
!
      open(unit=21,file=orbital_input_path(state_name),form="UNFORMATTED",         &
                                  status="OLD",position="asis")
!
!   Save NAK, NP and NH
!
      nak(:nwff) = nakff(:nwff)
      np(:nwff) = npff(:nwff)
      nh(:nwff) = nhff(:nwff)
!
!   Allocate storage to orbital arrays
!
      call alloc (pfff, max_radial_points,nwff, "PFFF", "LODRWFF")
      call alloc (qfff, max_radial_points,nwff, "QFFF", "LODRWF")
!
      con = nuclear_charge/speed_of_light_au
      con = con*con
!
      do j = 1, nwff
         pfff(:max_radial_points,j) = 0.0_real64
         qfff(:max_radial_points,j) = 0.0_real64
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
      read (21, iostat=io_status) radial_file_signature
      if (io_status/=0 .or. radial_file_signature/="G92RWF") then
         write (6, *) "This is not a Radial WaveFunction File;"
         close(21)
      end if

      if (ntest >= 100) then
         write (*, *) "******************"
         write (*, *) " Entering lodrwff"
         write (*, *) "******************"
      end if

    3 continue
      read (21, iostat=io_status) npy, naky, ey, my
      if (io_status == 0) then
         call alloc (pa, my, "PA", "LODRWFFF")
         call alloc (qa, my, "QA", "LODRWFFF")
         call alloc (ra, my, "RA", "LODRWFFF")
         read (21) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
         read (21) (ra(i),i=1,my)
!
!    Orbital order as defined in kapdata
!
         jj = 0
         do k = 1, kamax
            if (k > 1) jj = nshlff(k-1) + jj
            do j = 1, nshlff(k)
               kk = nshlpff(k,j)
               if (npy/=np(kk) .or. naky/=nak(kk)) cycle
               jjj = jj + j
               pzff(jjj) = pzy
               eff(jjj) = ey
               nakff(jjj) = nak(kk)
               npff(jjj) = np(kk)
               nhff(jjj) = nh(kk)
               gamaff(jjj) = gama(kk)
               call intrpqf (pa, qa, my, ra, jjj, dnorm)
               if (ntest >= 100) write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') npff(jjj), nhff(jjj), eff(jjj)&
                  , dnorm
               if (ntest > 1000) then
                  write (*, *) "PF              QF             RA"
                  do kkk = 1, mfff(jjj)
                     write (*, *) pfff(kkk,jjj), qfff(kkk,jjj), ra(kkk)
                  end do
               end if
               nwin = nwin + 1
            end do
         end do
         call dalloc (pa, "PA", "LODRWFF")
         call dalloc (qa, "QA", "LODRWFF")
         call dalloc (ra, "RA", "LODRWFF")

         goto 3
      end if
      if (ldbpr(3)) write (99, *) " orbitals renormalised;"
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nwff) then
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
         do j = 1, nwff
            write (*, '(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)') npff(j), nhff(j), eff(j), dnorm
         end do
         write (*, *)
         write (*, *) "*****************"
         write (*, *) " Leaving lodrwff"
         write (*, *) "*****************"
      end if

      return
!
      return
!
      end subroutine lodrwff
