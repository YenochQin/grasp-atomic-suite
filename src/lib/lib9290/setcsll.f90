!***********************************************************************
      subroutine setcsll(nunit, name, nblkin, nblock, ncfblk, ncftot, idblk)
!
!  Open, read name file to get nblock, ncfblk(), idblk(), ncftot
!
! Xinghong He 98-06-29
!
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:34   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use openfl_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: nunit
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
         write (6, *) name(1:len_trim(name)), ' does not exist'
         stop
      endif

! Open it

      call openfl (nunit, name, 'FORMATTED', 'OLD', ierr)
      if (ierr == 1) then
         write (6, *) 'Error when opening ', name(1:len_trim(name))
         stop
      endif

! Check the first record of the file; if not as expected, stop

      read (nunit, '(1A15)', iostat=ios) str
      if (ios/=0 .or. str/='Core subshells:') then
         write (6, *) 'Not a Configuration Symmetry List File;'
         close(nunit)
         stop
      endif

! Skip next 4 records

      do i = 1, 4
         read (nunit, *)
      end do

! Determine the number of blocks in this file

      nblock = 0
      ncsf = 0

      ios = 0
      do while(ios == 0)
         read (nunit, '(1A2)', iostat=ios) ch
         if (ch==' *' .or. ios/=0) then
            !.. a new block has been found
            nblock = nblock + 1
            write (6, *) 'Block ', nblock, ',  ncf = ', ncsf
            if (nblock > nblkin) then
               write (6, *) 'setcsll: Too many blocks(', nblock, ')'
               write (6, *) 'Maximum allowed is ', nblkin
               stop
            endif
            i = len_trim(line3)
            idblk(nblock) = line3(i-4:i)
            ncfblk(nblock) = ncsf
            ncsf = 0
            if (ios == 0) cycle
         else
            read (nunit, *)
            read (nunit, '(A)') line3
            ncsf = ncsf + 1
         endif
      end do

! Obtain ncftot

      ncftot = sum(ncfblk(:nblock))

      return
      end subroutine setcsll
