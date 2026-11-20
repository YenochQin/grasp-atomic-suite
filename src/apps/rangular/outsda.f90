!***********************************************************************
      subroutine outsda(print_details,nonzero_element_count, csf_count)
!***********************************************************************
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nonzero_element_count
      integer, intent(in) :: csf_count
      logical, intent(in) :: print_details
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=9), parameter :: routine_name ="outsdampi"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: nnonz_a
!-----------------------------------------------
      nnonz_a = nonzero_element_count

      write(6, *)" ... complete; density of non-zero elements of H(DC):", &
         nnonz_a, "/", (csf_count*(csf_count + 1))/2
!
!   Debug printout
!
      if (print_details) then
         write (99, *)
         write (99, *) "From ", routine_name, " :"
         write (99, '(" Number of nonzero elements in H(DC): ",1i4)') nnonz_a

         write (6, *) "This part not finished. See ", routine_name
         write (99, *) "This part not finished. See ", routine_name
      end if


      return
      end subroutine outsda
