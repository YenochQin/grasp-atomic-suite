

!***********************************************************************
!                                                                      *
      subroutine wrtrwf
!                                                                      *
!   Open, write a header and all subshell radial wavefunctions, and    *
!   close the  .rwf  file.                                             *
!                                                                      *
!   Call(s) to: [LIB92]: DALLOC, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!                                                                      *
!***********************************************************************
!
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use grid_c
      use orb_c
      use wave_c, only: orbital_origin_parameters, radial_large_components, radial_small_components, orbital_last_grid_index
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: newunit, ierr, j, mfj, i
      logical :: imopened
      character :: filnam*128
!-----------------------------------------------
!
!

      filnam = "rwfn.inp"

      do newunit = 23, 99                        ! 23 is a historical value
         inquire(unit=newunit, opened=imopened)
         if (imopened) cycle
         exit                                    ! should be the normal exit point
      end do

      if (newunit == 100) then
         write (stderr_unit, *) "All unit numbers from 23 to 99 are BUSY!"
         stop
      end if

      call openfl (newunit, filnam, "UNFORMATTED", "NEW", ierr)
      if (ierr == 1) then
         write (stderr_unit, *) 'Error when opening "', filnam(1:len_trim(filnam)), &
            '"'

         stop
      end if
!
!   Write the file header
!
      write (newunit) "G92RWF"
!
!   Write out the radial wavefunctions
!
      do j = 1, orbital_count
         mfj = orbital_last_grid_index(j)
         write (newunit) orbital_principal_n(j), orbital_kappa(j), orbital_energy_parameters(j), mfj
         write (newunit) orbital_origin_parameters(j), (radial_large_components(i,j),i=1,mfj), &
            & (radial_small_components(i,j),i=1,mfj)
         write (newunit) (radial_grid(i),i=1,mfj)
      end do

      close(newunit)
!
!   Deallocate the storage for the radial wavefunctions
!
      call dalloc (radial_large_components, "PF", "WRTRWF")
      call dalloc (radial_small_components, "QF", "WQRTRWF")
!
      return
      end subroutine wrtrwf
