!***********************************************************************
!                                                                      *
      subroutine setmixmpi(name, idblk)
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
      use mpi_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      use lodmixmpi_i
      use mpiu_i, only: stopmpi
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*) , intent(in) :: name
      character(len=8) :: idblk(*)
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=11), parameter :: form = 'UNFORMATTED'
      character(len=7), parameter  :: status = 'UNKNOWN'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: k
!-----------------------------------------------
!
      if (mpi_rank .eq. 0) then
         k = index (name,' ')
         call openfl (25, name(1:k-1)//'.cm', form, status, mpi_error_code)
         if (mpi_error_code .ne. 0) then
            call stopmpi ('setmix: Error when opening .cm file', mpi_rank)
         endif

         write (25) 'G92MIX'
      endif

      call lodmixmpi (idblk)

      return
      end subroutine setmixmpi
