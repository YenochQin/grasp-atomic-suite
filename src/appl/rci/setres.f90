!***********************************************************************
!                                                                      *
      subroutine setres(isofile, rwffile, idblk)
!                                                                      *
!   Open, check, load data from the  .res  file.                       *
!                                                                      *
!   Call(s) to: [LIB92]: GETYN, OPENFL.
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
!   Modified by Xinghong                  Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man, only: alloc, double
      use default_c, only: double, inpci, ndef, ndump
      use where_c, only: double, imcdf
      use hblock_c, only: csf_counts_by_block, double, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use iccu_c, only: iccut, iccutblk
      use iounit_c, only: stderr_unit, stdin_unit, stdout_unit
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use openfl_i, only: openfl
      use lodres_i, only: lodres
      use getcid_i, only: getcid
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: isofile*(*)
      character  :: rwffile*(*)
      character(len=8)  :: idblk(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=11), parameter :: form = "UNFORMATTED"
      character(len=7), parameter :: status = "UNKNOWN"
      character(len=6), parameter :: restitle = "R92RES"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ierr, ios
      logical :: found, restrt
      character :: r92res*6, defnam*11, idstring*3
!-----------------------------------------------
!
! Compose the "rci.res" file name
!
      defnam = "rci.res"
!
! Ask if this is a restart
!
      if (ndef /= 0) then
         write (stderr_unit, *) "Restarting RCI90 ?"
         restrt = getyn()
      else
         restrt = .false.
      end if
!      IF (RESTRT) THEN
!        WRITE(734,'(a)') 'y            ! Restarting RCI90 ?'
!      ELSE
!        WRITE(734,'(a)') 'n            ! Restarting RCI90 ?'
!      END IF
!
! Do some settings and checks
!
      if (restrt) then
!         ...Restart, make sure file exist
         inquire(file=defnam, exist=found)
         if (.not.found) stop "setres: .res does not exist"
      end if
!
! Open the .res file
!
      call openfl (imcdf, defnam, form, status, ierr)
      if (ierr /= 0) stop "setres: Error openning .res file"
!
! If restart, load the contents. Otherwise generate them via getcid
!
! But first of all, iccutblk() is needed in both cases
!
      call alloc (iccutblk, symmetry_block_count, "ICCUTBLK", "SETRES")

      if (restrt) then
!        ...Check the signature of the file
         read (imcdf, iostat=ios) r92res
         if (ios/=0 .or. r92res/=restitle) then
            close(imcdf)
            stop "setres: Not RCI92 .res file"
         end if

!         ...Read and check restart information
         call lodres

      else

!         ...Write the file header
!         ...Generate the first part of the .res file
         write (imcdf) restitle
         call getcid (isofile, rwffile, idblk)

      end if

      return
      end subroutine setres
