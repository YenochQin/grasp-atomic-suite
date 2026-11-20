!***********************************************************************
      subroutine outsda(lprint,nnonz, ncf)
!***********************************************************************
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nnonz
      integer, intent(in) :: ncf
      logical, intent(in) :: lprint
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=9), parameter :: myname ="outsdampi"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nnonz_a
!-----------------------------------------------
      nnonz_a = nnonz

      write(6, *)" ... complete; density of non-zero elements of H(DC):", &
         nnonz_a, "/", (ncf*(ncf + 1))/2
!
!   Debug printout
!
      if (lprint) then
         write (99, *)
         write (99, *) "From ", myname, " :"
         write (99, '(" Number of nonzero elements in H(DC): ",1i4)') nnonz_a

         write (6, *) "This part not finished. See ", myname
         write (99, *) "This part not finished. See ", myname
      end if


      return
      end subroutine outsda
