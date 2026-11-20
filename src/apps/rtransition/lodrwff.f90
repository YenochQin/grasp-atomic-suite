!***********************************************************************
!                                                                      *
      subroutine lodrwff(state_name)
      use, intrinsic :: iso_fortran_env, only: output_unit
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Written by Per Jonsson                              June 1996      *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only: max_radial_points
      use array_memory_management, only: alloc, dalloc
      use atomic_calculation_state, only: max_radial_points, nuclear_charge, speed_of_light_au
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
      integer :: wavefunction_input_unit
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: state_name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, io_status, npyff, nakyff, myff, jj
      real(real64) :: con, fkk, eyff, pzy
      real(real64), pointer :: pa(:), qa(:), ra(:)
      character :: radial_file_signature*6
!-----------------------------------------------
!
!   Common relevant for the final state
!
!
!   Write entry message
!
      write (output_unit, *) "Loading Radial WaveFunction File for final state..."
!
!   Open the radial wave function file
!
      j = index(state_name," ")
      open(newunit=wavefunction_input_unit, file=state_name(1:j-1)//".bw", form="UNFORMATTED", status="OLD", &
         position="asis")
!
!   Allocate storage to orbital arrays
!
      call alloc (pfff, max_radial_points, nwff, "PFFF", "LODRWFF")
      call alloc (qfff, max_radial_points, nwff, "QFFF", "LODRWFF")
!
!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      con = nuclear_charge/speed_of_light_au
      con = con*con
!
      do j = 1, nwff
         pfff(:max_radial_points,j) = 0.0_real64
         qfff(:max_radial_points,j) = 0.0_real64
!
         eff(j) = -1.0_real64
!
         k = abs(nakff(j))
         if (nuclear_model_parameter_count > 0) then
            gamaff(j) = dble(k)
         else if (nuclear_model_parameter_count == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               gamaff(j) = sqrt(fkk - con)
            else
               write (output_unit, *) "LODRWF: Imaginary gamma parameter"
               write (output_unit, *) " for ", npff(j), nhff(j), " orbital; the"
               write (output_unit, *) " point model for the nucleus"
               write (output_unit, *) " is inappropriate for Z > ", speed_of_light_au, "."
               stop
            end if
         end if
!
      end do
!
!   Read orbital information from Read Orbitals File;
!
      nwin = 0
      read (wavefunction_input_unit, iostat=io_status) radial_file_signature
      if (io_status/=0 .or. radial_file_signature/="G92RWF") then
         write (output_unit, *) "This is not a Radial WaveFunction File;"
         close(wavefunction_input_unit)
      end if

    3 continue
      read (wavefunction_input_unit, iostat=io_status) npyff, nakyff, eyff, myff
      if (io_status == 0) then
         call alloc (pa, myff, "PA", "LODRWFF" )
         call alloc (qa, myff, "QA", "LODRWFF" )
         call alloc (ra, myff, "RA", "LODRWFF" )
         read (wavefunction_input_unit) pzy, (pa(i),i=1,myff), (qa(i),i=1,myff)
         read (wavefunction_input_unit) (ra(i),i=1,myff)

         do j = 1, nwff
            if (.not.(eff(j)<0.0_real64 .and. npyff==npff(j) .and. nakyff==nakff(j)&
               )) cycle
            pzff(j) = pzy
            eff(j) = eyff
            mfff(j) = myff
            do jj = 1, mfff(j)
               pfff(jj,j) = pa(jj)
               qfff(jj,j) = qa(jj)
            end do
            nwin = nwin + 1
         end do
         call dalloc (pa, "PA", "LODRWFF")
         call dalloc (qa, "QA", "LODRWFF")
         call dalloc (ra, "RA", "LODRWFF")
         goto 3
      end if
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nwff) then
         write (output_unit, *) "LODRWF: All required orbitals not"
         write (output_unit, *) " found."
         stop
      end if
!
      write (output_unit, *) " ... load complete;"
!
      close(wavefunction_input_unit)

      return
      end subroutine lodrwff
