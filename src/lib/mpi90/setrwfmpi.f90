module setrwfmpi_i
    implicit none (type, external)
    private
    public :: setrwfmpi
contains
!***********************************************************************
      subroutine setrwfmpi (name)

!   Open, check, load data from and close the  .rwf  file.             *
!
!   Used by rcimpivu, rscfmpi, rcimpi
!                                                                      *
!   Call(s) to: [LIB92]: LODRWFmpi, OPENFL.                            *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   MPI version by Xinghong He            Last revision: 06 Aug 1998   *
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_c, only: mpi_rank
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i, only: openfl
      use lodrwfmpi_i, only: lodrwfmpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = *) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ios, ierror
      character (len = 6) :: g92rwf
!-----------------------------------------------
      if (mpi_rank == 0) then
         call openfl (23, name, "UNFORMATTED", "OLD", ierror)
         if (ierror == 1) then
            write (stderr_unit,*) "Error opening", name(1:len_trim (name))
            stop
         end if

!   Check the file; if not as expected, stop.

         read (23,iostat = ios) g92rwf
         if ((ios /= 0) .or. (g92rwf /= "G92RWF")) then
            write (stderr_unit,*) "This is not a Radial WaveFunction File;"
            close (23)
            stop
         end if
      end if

!   Attempt to load the radial wavefunctions

      call lodrwfmpi (ierror)

      if (ierror /= 0) then
         if (mpi_rank == 0) then
            write (stderr_unit,*) "Radial wavefunctions defined in CSL file" &
                          , " not found in Radial WaveFunction File"
            close (23)
         end if
         stop
      end if

      if (mpi_rank == 0) close (23)

      return
      end subroutine setrwfmpi
end module setrwfmpi_i
