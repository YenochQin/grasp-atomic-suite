module setcsll_i
    implicit none (type, external)
    private
    public :: setcsll
contains
!***********************************************************************
      subroutine setcsll(nunit, name, nblkin, nblock, ncfblk, ncftot, idblk)
      use, intrinsic :: iso_fortran_env, only: output_unit
!
!  Open, read name file to get nblock, ncfblk(), idblk(), ncftot
!
! Xinghong He 98-06-29
!
!***********************************************************************
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use file_open, only: openfl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: nunit
      integer, intent(in) :: nblkin
      integer, intent(out) :: nblock
      integer, intent(out) :: ncftot
      character (len = *), intent(inout) :: name
      integer, intent(inout) :: ncfblk(*)
      character (len = 8), intent(out) :: idblk(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, ncsf, ios, ierr
      logical :: found
      character :: str*15, ch*2, line3*200
!-----------------------------------------------
! Locals

! Look for  <name>

      inquire(file=name, exist=found)
      if (.not.found) then
         write (output_unit, *) name(1:len_trim(name)), " does not exist"
         stop
      end if

! Open it

      call openfl (nunit, name, "FORMATTED", "OLD", ierr)
      if (ierr == 1) then
         write (output_unit, *) "Error when opening ", name(1:len_trim(name))
         stop
      end if

! Check the first record of the file; if not as expected, stop

      read (nunit, "(1A15)", iostat=ios) str
      if (ios/=0 .or. str/="Core subshells:") then
         write (output_unit, *) "Not a Configuration Symmetry List File;"
         close(nunit)
         stop
      end if

! Skip next 4 records

      do i = 1, 4
         read (nunit, *)
      end do

! Determine the number of blocks in this file

      nblock = 0
      ncsf = 0

      ios = 0
      do while(ios == 0)
         read (nunit, "(1A2)", iostat=ios) ch
         if (ch==" *" .or. ios/=0) then
            !.. a new block has been found
            nblock = nblock + 1
            write (output_unit, *) "Block ", nblock, ",  ncf = ", ncsf
            if (nblock > nblkin) then
               write (output_unit, *) "setcsll: Too many blocks(", nblock, ")"
               write (output_unit, *) "Maximum allowed is ", nblkin
               stop
            end if
            i = len_trim(line3)
            idblk(nblock) = line3(i-4:i)
            ncfblk(nblock) = ncsf
            ncsf = 0
            if (ios == 0) cycle
         else
            read (nunit, *)
            read (nunit, "(A)") line3
            ncsf = ncsf + 1
         end if
      end do

! Obtain ncftot

      ncftot = sum(ncfblk(:nblock))

      return
      end subroutine setcsll
end module setcsll_i
