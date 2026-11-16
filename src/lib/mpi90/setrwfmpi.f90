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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use mpi_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodrwfmpi_i
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
      if (mpi_rank .eq. 0) then
         call openfl (23, name, 'UNFORMATTED', 'OLD', ierror)
         if (ierror .eq. 1) then
            write (stderr_unit,*) 'Error opening', name(1:len_trim (name))
            stop
         endif

!   Check the file; if not as expected, stop.

         read (23,iostat = ios) g92rwf
         if ((ios .ne. 0) .or. (g92rwf .ne. 'G92RWF')) then
            write (stderr_unit,*) 'This is not a Radial WaveFunction File;'
            close (23)
            stop
         endif
      endif

!   Attempt to load the radial wavefunctions

      call lodrwfmpi (ierror)

      if (ierror .ne. 0) then
         if (mpi_rank .eq. 0) then
            write (stderr_unit,*) 'Radial wavefunctions defined in CSL file' &
                          , ' not found in Radial WaveFunction File'
            close (23)
         endif
         stop
      endif

      if (mpi_rank .eq. 0) close (23)

      return
      end
