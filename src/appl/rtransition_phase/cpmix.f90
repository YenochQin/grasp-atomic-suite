!************************************************************************
      subroutine cpmix(name, inpci)
!************************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: inpci
      character , intent(inout) :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ios, nelec, ncftot, nw, nvectot, nvecsiz, nblock, iblk, ib, &
         ncf, nvec, iatjp, iaspa, i, j
      real(real64), pointer :: eval(:), evec(:)
      integer, pointer :: ivec(:)
      real(real64) :: eav
      character :: g92mix*6
!-----------------------------------------------


      name(2) = trim(name(1))//'_CP'
      if (inpci == 0) then
         open(unit=78, file=trim(name(2))//'.cbm', form='UNFORMATTED', status=&
            'UNKNOWN', position='asis')
      else
         open(unit=78, file=trim(name(2))//'.bm', form='UNFORMATTED', status=&
            'UNKNOWN', position='asis')
      endif
      if (inpci == 0) then
         open(unit=68, file=trim(name(1))//'.cbm', form='UNFORMATTED', status=&
            'OLD', position='asis')
      else
         open(unit=68, file=trim(name(1))//'.bm', form='UNFORMATTED', status=&
            'OLD', position='asis')
      endif
      read (68, iostat=ios) g92mix
      if (ios/=0 .or. g92mix/='G92MIX') then
         write (*, *) 'Not a GRASP mixing file'
         stop
      endif
      write (78) g92mix
      read (68) nelec, ncftot, nw, nvectot, nvecsiz, nblock
      write (78) nelec, ncftot, nw, nvectot, nvecsiz, nblock

      do iblk = 1, nblock
         read (68) ib, ncf, nvec, iatjp, iaspa
         write (78) ib, ncf, nvec, iatjp, iaspa
         write (*,*) ib, ncf, nvec, iatjp, iaspa
         call alloc (eval, nvec, 'EVAL', 'CPMIX')
         call alloc (evec, ncf*nvec, 'EVEC', 'CPMIX')
         call alloc (ivec, nvec, 'IVEC', 'CPMIX')
         read (68) (ivec(i),i=1,nvec)
         read (68) eav, (eval(i),i=1,nvec)
         read (68) ((evec(i + (j - 1)*ncf),i=1,ncf),j=1,nvec)

         write (78) (ivec(i),i=1,nvec)
         write (78) eav, (eval(i),i=1,nvec)
         write (78) ((evec(i + (j - 1)*ncf),i=1,ncf),j=1,nvec)
         call dalloc (eval, 'EVAL', 'CPMIX')
         call dalloc (evec, 'EVEC', 'CPMIX')
         call dalloc (ivec, 'IVEC', 'CPMIX')
      end do
      name(2) = name(1)
      close(68)
      close(78)
      return
      end subroutine cpmix
