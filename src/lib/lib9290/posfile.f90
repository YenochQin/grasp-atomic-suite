module posfile_i
    implicit none (type, external)
    private
    public :: posfile
contains
!***********************************************************************
      subroutine posfile(mode, nunit, nrec)
!    Position the file
!***********************************************************************
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
end module posfile_i
