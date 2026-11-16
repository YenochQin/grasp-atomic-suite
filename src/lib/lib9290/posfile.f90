!***********************************************************************
      subroutine posfile(mode, nunit, nrec)
!    Position the file
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:12   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!----------------------------------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: mode
      integer, intent(in) :: nunit
      integer, intent(in) :: nrec
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: forwd = 0
      integer, parameter :: backwd = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
!

      select case (mode)

      case (forwd)
         rewind (nunit)
         do i = 1, nrec
            read (nunit)
         end do

      case (backwd)
         do i = 1, nrec
            backspace (nunit)
         end do

      end select

      return
      end subroutine posfile
