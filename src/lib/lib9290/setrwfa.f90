!***********************************************************************
!                                                                      *
      subroutine setrwfa(name)
!                                                                      *
!   Open, check, load data from and close the  .rwf  file.             *
!                                                                      *
!   Call(s) to: [LIB92]: LODRWF, OPENFL.                               *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified    by Xinghong He            Last revision: 09 Jul 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodrwf_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = *), intent(in) :: name
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr, ios
      character :: g92rwf*6


      call openfl (23, name, 'UNFORMATTED', 'OLD', ierr)
      if (ierr == 1) then
         write (stderr_unit, *) 'Error when opening', name(1:len_trim(name))
         stop
      endif
!
!   Check the file; if not as expected, stop.
!
      read (23, iostat=ios) g92rwf
      if (ios/=0 .or. g92rwf/='G92RWF') then
         write (stderr_unit, *) 'This is not a Radial WaveFunction File;'
         close(23)
         stop
      endif
!
!   Attempt to load the radial wavefunctions; if this fails, stop
!
      call lodrwf (ierr)

      if (ierr /= 0) then
         write (stderr_unit, *) 'Radial wavefunctions defined in CSL file', &
            ' not found in Radial WaveFunction File'
         close(23)
         stop
      endif

      close(23)

      return
      end subroutine setrwfa
