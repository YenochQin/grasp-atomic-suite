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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use memory_man
      use default_c
      use where_c
      use hblock_c
      use iccu_c
      use iounit_c
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use openfl_i
      use lodres_i
      use getcid_i
      use mpiu_i, only: stopmpi
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: isofile*(*)
      character  :: rwffile*(*)
      character(len=8), dimension(*)  :: idblk
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character*11, parameter :: form = 'UNFORMATTED'
      character*7, parameter :: status = 'UNKNOWN'
      character*6, parameter :: restitle = 'R92RES'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ios
      logical :: found, restrt
      character :: r92res*6, defnam*11, idstring*3
!-----------------------------------------------
!
! Compose the "rciXXX.res" file names (for each node)
!
      write (idstring, '(I3.3)') mpi_rank
      defnam = 'rci' // idstring // '.res'
!
! Ask if this is a restart
!
      if (mpi_rank .eq. 0) then
         if (ndef /= 0) then
            write (stderr_unit, *) 'Restarting RCI90 ?'
            restrt = getyn()
         else
            restrt = .false.
         endif
!         IF (RESTRT) THEN
!           WRITE(734,'(a)') 'y            ! Restarting RCI90 ?'
!         ELSE
!           WRITE(734,'(a)') 'n            ! Restarting RCI90 ?'
!         END IF
      end if
      call mpi_bcast (restrt, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
!
! Do some settings and checks
!
      if (restrt) then
!         ...Restart, make sure file exist
         inquire(file=defnam, exist=found)
         if (.not. found) then
            call stopmpi ('setres: .res not exist', mpi_rank)
         endif
      endif
!
! Open the .res file
!
      call openfl (imcdf, defnam, form, status, mpi_error_code)
      if (mpi_error_code .ne. 0) then
         call stopmpi ('setres: Error openning .res file', mpi_rank)
      endif
!
! If restart, load the contents. Otherwise generate them via getcid
!
! But first of all, iccutblk() is needed in both cases
!
      call alloc (iccutblk, symmetry_block_count, 'ICCUTBLK', 'SETRES')

      if (restrt) then
!        ...Check the signature of the file
         read (imcdf, iostat=ios) r92res
         if (ios/=0 .or. r92res/=restitle) then
            close(imcdf)
            call stopmpi ('setres: Not RCI92 .res file', mpi_rank)
         endif

!         ...Read and check restart information
         call lodres

      else

!         ...Write the file header
!         ...Generate the first part of the .res file
         write (imcdf) restitle
         call getcid (isofile, rwffile, idblk)

      endif

      return
      end subroutine setres
