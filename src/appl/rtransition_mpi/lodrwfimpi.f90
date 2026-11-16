!***********************************************************************
!                                                                      *
      subroutine lodrwfi(name)
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
      use parameter_def,   only: max_radial_points
      use memory_man
      use def_c
      use npar_c, only: double, nuclear_model_parameter_count, nuclear_model_parameters
      use biorb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: name*128
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npyii, nakyii, myii, jj
      real(real64) :: con, fkk, eyii, pzy
      real(real64), pointer :: pa(:), qa(:), ra(:)
      character :: g92rwf*6
!-----------------------------------------------
!
!   Common relevant for the initial state
!
!   Write entry message
!
!      WRITE (6, *) 'Loading Radial WaveFunction File for initial state...'
!
!   Open the radial wave function file
!
      j = index(name," ")
      open(unit=69, file=name(1:j-1)//".bw", form="UNFORMATTED", status="OLD", &
         position="asis")
!
!   Allocate storage to orbital arrays
!
      call alloc (pfii, max_radial_points, nwii, "PFII", "LODRWFI")
      call alloc (qfii, max_radial_points, nwii, "QFII", "LODRWFI")
!
      con = nuclear_charge/speed_of_light_au
      con = con*con
!
      write (*, *) "NWII", nwii
      do j = 1, nwii
         write (*, *) nakii(j), npii(j), nhii(j)
         pfii(:max_radial_points,j) = 0.0_real64
         qfii(:max_radial_points,j) = 0.0_real64
!
         eii(j) = -1.0_real64
!
         k = abs(nakii(j))
         if (nuclear_model_parameter_count > 0) then
            gamaii(j) = dble(k)
         else if (nuclear_model_parameter_count == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               gamaii(j) = sqrt(fkk - con)
            else
               write (6, *) "LODRWF: Imaginary gamma parameter"
               write (6, *) " for ", npii(j), nhii(j), " orbital; the"
               write (6, *) " point model for the nucleus"
               write (6, *) " is inappropriate for Z > ", speed_of_light_au, "."
               stop
            end if
         end if
!
      end do
!
!   Read orbital information from Read Orbitals File;
!
      nwin = 0
      read (69, iostat=ios) g92rwf
      if (ios/=0 .or. g92rwf/="G92RWF") then
         write (6, *) "This is not a Radial WaveFunction File;"
         close(69)
      end if

    3 continue
      read (69, iostat=ios) npyii, nakyii, eyii, myii
      if (ios == 0) then
         call alloc (pa, myii, "PA", "LODRWFI")
         call alloc (qa, myii, "QA", "LODRWFI")
         call alloc (ra, myii, "RA", "LODRWFI")
         read (69) pzy, (pa(i),i=1,myii), (qa(i),i=1,myii)
         read (69) (ra(i),i=1,myii)

         do j = 1, nwii
            if (.not.(eii(j)<0.0_real64 .and. npyii==npii(j) .and. nakyii==nakii(j)&
               )) cycle
            pzii(j) = pzy
            eii(j) = eyii
            mfii(j) = myii
            do jj = 1, mfii(j)
               pfii(jj,j) = pa(jj)
               qfii(jj,j) = qa(jj)
            end do
            nwin = nwin + 1
         end do
         call dalloc (pa, "PA", "LODRWFI")
         call dalloc (qa, "QA", "LODRWFI")
         call dalloc (ra, "RA", "LODRWFI")
         goto 3
      end if
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nwii) then
         write (6, *) "LODRWF: All required orbitals not"
         write (6, *) " found."
         stop
      end if
!
!      WRITE (6, *) ' ... load complete;'
!
      close(69)

      return
      end subroutine lodrwfi
