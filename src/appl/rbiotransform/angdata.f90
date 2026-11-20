!***********************************************************************
!                                                                      *
      subroutine angdata(name, avail, kamax)
!                                                                      *
!   Checks if the angular file name.T is available and appropriate     *
!                                                                      *
!   Written by Per Jonsson                      6 March 1997           *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis,           only: csf_count, orbital_count
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)      :: kamax
      logical, intent(out)     :: avail
      character, intent(inout) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, nf, ncfd, nwd, kamaxd
      logical :: found
!-----------------------------------------------
!
      j = index(name," ")
      inquire(file=name(1:j-1)//".TB", exist=found)
      if (.not.found) then
         write (6, *) " Angular file not available"
         avail = .false.
         return
      else
!
!  Open the file and check if it is appropriate for the present case
!
         nf = 200
         open(unit=nf, file=name(1:j-1)//".TB", status="OLD", form=&
            "UNFORMATTED", position="asis")
         rewind (nf)
         read (nf) ncfd, nwd, kamaxd
         if (.not.(ncfd==csf_count .and. nwd==orbital_count .and. kamaxd==kamax)) then
            write (6, *) " Angular file not appropriate"
            avail = .false.
            return
         else
            write (6, *) " Angular data read from file"
            avail = .true.
         end if
      end if
      return
      end subroutine angdata
