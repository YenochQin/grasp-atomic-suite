!***********************************************************************
!                                                                      *
      subroutine setres(isotope_file_path, radial_file_path, block_labels)
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
      use array_memory_management, only: alloc
      use calculation_input_defaults, only: inpci, ndef, ndump
      use where_c, only: imcdf
      use hamiltonian_symmetry_blocks, only: csf_counts_by_block, eigenvector_counts_by_block, ncmaxblk, symmetry_block_count
      use csf_block_cutoffs, only: iccut, iccutblk
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use mpi_state, only: mpi_bcast, mpi_comm_world, mpi_error_code, mpi_logical, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use file_open, only: openfl
      use lodres_i, only: lodres
      use getcid_i, only: getcid
      use mpi_utilities, only: stopmpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character  :: isotope_file_path*(*)
      character  :: radial_file_path*(*)
      character(len=8)  :: block_labels(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=11), parameter :: form = "UNFORMATTED"
      character(len=7), parameter :: status = "UNKNOWN"
      character(len=6), parameter :: restitle = "R92RES"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: io_status
      logical :: found, is_restart
      character :: r92res*6, file_path*11, block_label_text*3
!-----------------------------------------------
!
! Compose the "rciXXX.res" file names (for each node)
!
      write (block_label_text, "(I3.3)") mpi_rank
      file_path = "rci" // block_label_text // ".res"
!
! Ask if this is a restart
!
      if (mpi_rank == 0) then
         if (ndef /= 0) then
            write (stderr_unit, *) "Restarting RCI90 ?"
            is_restart = getyn()
         else
            is_restart = .false.
         end if
!         IF (RESTRT) THEN
!           WRITE(734,'(a)') 'y            ! Restarting RCI90 ?'
!         ELSE
!           WRITE(734,'(a)') 'n            ! Restarting RCI90 ?'
!         END IF
      end if
      call mpi_bcast (is_restart, 1, mpi_logical, 0, mpi_comm_world, mpi_error_code)
!
! Do some settings and checks
!
      if (is_restart) then
!         ...Restart, make sure file exist
         inquire(file=file_path, exist=found)
         if (.not. found) then
            call stopmpi ("setres: .res not exist", mpi_rank)
         end if
      end if
!
! Open the .res file
!
      call openfl (imcdf, file_path, form, status, mpi_error_code)
      if (mpi_error_code /= 0) then
         call stopmpi ("setres: Error openning .res file", mpi_rank)
      end if
!
! If restart, load the contents. Otherwise generate them via getcid
!
! But first of all, iccutblk() is needed in both cases
!
      call alloc (iccutblk, symmetry_block_count, "ICCUTBLK", "SETRES")

      if (is_restart) then
!        ...Check the signature of the file
         read (imcdf, iostat=io_status) r92res
         if (io_status/=0 .or. r92res/=restitle) then
            close(imcdf)
            call stopmpi ("setres: Not RCI92 .res file", mpi_rank)
         end if

!         ...Read and check restart information
         call lodres

      else

!         ...Write the file header
!         ...Generate the first part of the .res file
         write (imcdf) restitle
         call getcid (isotope_file_path, radial_file_path, block_labels)

      end if

      return
      end subroutine setres
