!***********************************************************************
      subroutine outsda(lprint,nnonz, ncf)
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  11:01:42   1/ 5/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
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
      character(len=9), parameter :: myname ='outsdampi'
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nnonz_a
!-----------------------------------------------
      nnonz_a = nnonz

      write(6, *)' ... complete; density of non-zero elements of H(DC):', &
         nnonz_a, '/', (ncf*(ncf + 1))/2
!
!   Debug printout
!
      if (lprint) then
         write (99, *)
         write (99, *) 'From ', myname, ' :'
         write (99, 301) nnonz_a

         write (6, *) 'This part not finished. See ', myname
         write (99, *) 'This part not finished. See ', myname
      endif

  301 format(' Number of nonzero elements in H(DC): ',1i4)
  302 format(' Column ',1i2,', row ',1i2,', sparse matrix index ',1i4)

      return
      end subroutine outsda
