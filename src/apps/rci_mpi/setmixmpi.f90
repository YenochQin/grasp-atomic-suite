!***********************************************************************
!                                                                      *
      subroutine setmixmpi(name, idblk)
      use calculation_state_files, only: ci_mixing_suffix
!                                                                      *
!   Opens the  .mix  file on stream 25; writes a header to this file;  *
!   calls LODMIX to interactively determine the eigenpairs required.   *
!                                                                      *
!   Call(s) to: [LIB92]: OPENFL.                                       *
!               [RCI92]: LODMIX.                                       *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
!   Modified by Xinghong He               Last revision: 23 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      use mpi_state, only: mpi_error_code, mpi_rank
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      use lodmixmpi_i, only: lodmixmpi
      use mpi_utilities, only: stopmpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: name
      character(len=8) :: idblk(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=11), parameter :: form = "UNFORMATTED"
      character(len=7), parameter  :: status = "UNKNOWN"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
      if (mpi_rank == 0) then
         call openfl (25, trim(name)//ci_mixing_suffix, form, status, mpi_error_code)
         if (mpi_error_code /= 0) then
            call stopmpi ("setmix: Error when opening .cimix file", mpi_rank)
         end if

         write (25) "G92MIX"
      end if

      call lodmixmpi (idblk)

      return
      end subroutine setmixmpi
