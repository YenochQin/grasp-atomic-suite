!***********************************************************************
!                                                                      *
      subroutine lodrwff(name)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Written by Per Jonsson                              June 1996      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  07:25:11   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use memory_man
      use def_c
      use npar_c
      use biorb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character , intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npyff, nakyff, myff, jj
      real(double) :: con, fkk, eyff, pzy
      real(double), dimension(:), pointer :: pa, qa, ra
      character :: g92rwf*6
!-----------------------------------------------
!
!   Common relevant for the final state
!
!
!   Write entry message
!
      write (6, *) 'Loading Radial WaveFunction File for final state...'
!
!   Open the radial wave function file
!
      j = index(name,' ')
      open(unit=69, file=name(1:j-1)//'.bw', form='UNFORMATTED', status='OLD', &
         position='asis')
!
!   Allocate storage to orbital arrays
!
      call alloc (pfff, nnnp, nwff, 'PFFF', 'LODRWFF')
      call alloc (qfff, nnnp, nwff, 'QFFF', 'LODRWFF')
!
!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      con = z/c
      con = con*con
!
      do j = 1, nwff
         pfff(:nnnp,j) = 0.0d00
         qfff(:nnnp,j) = 0.0d00
!
         eff(j) = -1.0d00
!
         k = abs(nakff(j))
         if (nparm > 0) then
            gamaff(j) = dble(k)
         else if (nparm == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               gamaff(j) = sqrt(fkk - con)
            else
               write (6, *) 'LODRWF: Imaginary gamma parameter'
               write (6, *) ' for ', npff(j), nhff(j), ' orbital; the'
               write (6, *) ' point model for the nucleus'
               write (6, *) ' is inappropriate for Z > ', c, '.'
               stop
            endif
         endif
!
      end do
!
!   Read orbital information from Read Orbitals File;
!
      nwin = 0
      read (69, iostat=ios) g92rwf
      if (ios/=0 .or. g92rwf/='G92RWF') then
         write (6, *) 'This is not a Radial WaveFunction File;'
         close(69)
      endif

    3 continue
      read (69, iostat=ios) npyff, nakyff, eyff, myff
      if (ios == 0) then
         call alloc (pa, myff, 'PA', 'LODRWFF' )
         call alloc (qa, myff, 'QA', 'LODRWFF' )
         call alloc (ra, myff, 'RA', 'LODRWFF' )
         read (69) pzy, (pa(i),i=1,myff), (qa(i),i=1,myff)
         read (69) (ra(i),i=1,myff)

         do j = 1, nwff
            if (.not.(eff(j)<0.0d00 .and. npyff==npff(j) .and. nakyff==nakff(j)&
               )) cycle
            pzff(j) = pzy
            eff(j) = eyff
            mfff(j) = myff
            do jj = 1, mfff(j)
               pfff(jj,j) = pa(jj)
               qfff(jj,j) = qa(jj)
            end do
            nwin = nwin + 1
         end do
         call dalloc (pa, 'PA', 'LODRWFF')
         call dalloc (qa, 'QA', 'LODRWFF')
         call dalloc (ra, 'RA', 'LODRWFF')
         go to 3
      endif
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nwff) then
         write (6, *) 'LODRWF: All required orbitals not'
         write (6, *) ' found.'
         stop
      endif
!
      write (6, *) ' ... load complete;'
!
      close(69)

      return
      end subroutine lodrwff
