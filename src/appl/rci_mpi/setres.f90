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
      use mpi_c, only: double, mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i, only: getyn
      use openfl_i, only: openfl
      use lodres_i, only: lodres
      use getcid_i, only: getcid
      use mpiu_i, only: stopmpi
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
      integer :: ios
      logical :: found, restrt
      character :: r92res*6, defnam*11, idstring*3
!-----------------------------------------------
!
! Compose the "rciXXX.res" file names (for each node)
!
      write (idstring, "(I3.3)") mpi_rank
      defnam = "rci" // idstring // ".res"
!
! Ask if this is a restart
!
      if (mpi_rank == 0) then
         if (ndef /= 0) then
            write (stderr_unit, *) "Restarting RCI90 ?"
            restrt = getyn()
         else
            restrt = .false.
         end if
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
            call stopmpi ("setres: .res not exist", mpi_rank)
         end if
      end if
!
! Open the .res file
!
      call openfl (imcdf, defnam, form, status, mpi_error_code)
      if (mpi_error_code /= 0) then
         call stopmpi ("setres: Error openning .res file", mpi_rank)
      end if
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
            call stopmpi ("setres: Not RCI92 .res file", mpi_rank)
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
