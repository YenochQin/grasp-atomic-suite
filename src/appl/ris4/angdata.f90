!***********************************************************************
!                                                                      *
      subroutine angdata(name,avail,whichone)
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
      use vast_kind_param,  only: double
      use orb_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character (len = 24), intent(in) :: name
      logical, intent(out) :: avail
      integer, intent(in)  :: whichone
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      logical :: found
      integer :: j, nf, istat, iword
!-----------------------------------------------
!
      j = index(name,' ')
      if (whichone.eq.1) then
          nf = 50
          inquire (file = name(1:j-1)//'.IOB', exist = found)
      elseif (whichone.eq.2) then
          nf = 51
          inquire (file = name(1:j-1)//'.ITB', exist = found)
      else
         print*,'This message should never appear'
      endif

      if (.not.found) then
      if (whichone.eq.1) print *, ' One-body angular file not available'
      if (whichone.eq.2) print *, ' Two-body angular file not available'
        avail = .false.
        return
      else
        if (whichone.eq.1) then
          open(unit=nf,file = name(1:j-1)//'.IOB',status='UNKNOWN'     &
            ,position='APPEND',form='UNFORMATTED')
        else
          open(unit=nf,file = name(1:j-1)//'.ITB',status='UNKNOWN'     &
            ,position='APPEND',form='UNFORMATTED')
        endif
        backspace (unit = nf,iostat=istat)
        read(nf) iword
        rewind(unit = nf)
        if (iword.eq.-1) then
          if (whichone.eq.1) print *, ' One-body angular file available'
          if (whichone.eq.2) print *, ' Two-body angular file available'
          avail = .true.
        else
          avail = .false.
          if (whichone.eq.1) print *,' One-body angular file incomplete'
          if (whichone.eq.2) print *,' Two-body angular file incomplete'
        endif
      endif
      return
      end subroutine angdata
