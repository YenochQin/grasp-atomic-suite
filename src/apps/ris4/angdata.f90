!***********************************************************************
!                                                                      *
      subroutine angdata(name,avail,whichone,status)
!                                                                      *
!   Checks if the angular files name.IOB (one-body) or                 *
!   name.ITB (Two-body) are available                                  *
!                                                                      *
!                                                                      *
!   Modified by C. Naz\'e  Oct. 2011                                   *
!                                                                      *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: iostat_end
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = 24), intent(in) :: name
      logical, intent(out) :: avail
      integer, intent(in)  :: whichone
      integer, intent(out) :: status
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: found
      integer :: nf, istat, iword, close_status
      character(len=28) :: filename
      character :: probe
!-----------------------------------------------
!
      avail = .false.
      status = 0
      if (whichone==1) then
          nf = 50
          filename = trim(name)//".IOB"
      else if (whichone==2) then
          nf = 51
          filename = trim(name)//".ITB"
      else
         status = 1
         return
      end if
      inquire(file=trim(filename), exist=found, iostat=status)
      if (status /= 0) return
      if (.not.found) then
      if (whichone==1) print *, " One-body angular file not available"
      if (whichone==2) print *, " Two-body angular file not available"
        avail = .false.
        return
      end if
      ! Missing/incomplete optional caches may be recomputed. Access and
      ! positioning failures on an otherwise usable cache are fatal to callers.
      ! Check raw readability separately from sequential-record decoding so an
      ! OS read error is not mistaken for an incomplete record (e.g. a directory).
      open(unit=nf,file=trim(filename),status="OLD",action="READ", &
           access="STREAM",form="UNFORMATTED",iostat=status)
      if (status /= 0) return
      read(nf,iostat=status) probe
      if (status == iostat_end) status = 0
      close(unit=nf,iostat=close_status)
      if (close_status /= 0 .and. status == 0) status = close_status
      if (status /= 0) return
      open(unit=nf,file=trim(filename),status="OLD",action="READ", &
           position="APPEND",form="UNFORMATTED",iostat=status)
      if (status /= 0) return
      backspace(unit=nf,iostat=istat)
      if (istat == 0) then
        read(nf,iostat=istat) iword
        if (istat == 0) avail = iword == -1
      end if
      if (avail) then
        rewind(unit=nf,iostat=status)
        if (status == 0) then
          if (whichone==1) print *, " One-body angular file available"
          if (whichone==2) print *, " Two-body angular file available"
          return
        end if
      end if
      avail = .false.
      close(unit=nf,iostat=close_status)
      if (close_status /= 0 .and. status == 0) status = close_status
      if (status /= 0) return
      print *, " Angular file incomplete: ", trim(filename), "; iostat=", istat
      return
      end subroutine angdata
