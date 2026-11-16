!***********************************************************************
!                                                                      *
      subroutine lodrwfi(name, ntestg)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Per Jonsson                              June 1996      *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:29:50   1/ 6/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,   only: max_orbitals, max_radial_points
      use memory_man
      use biorb_c
      use def_c,           only: nuclear_charge, speed_of_light_au
      use debug_c
      use grid_c
      use npar_c
      use sbdat_c,         only: kamax, nshlii, nshlpii
      use wave_c,          only: pfii, qfii, pzii, mfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpqi_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer , intent(in) :: ntestg
      character , intent(in) :: name*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer  :: nak(max_orbitals), np(max_orbitals)
      integer :: ntestl, ntest, j, k, i, nwin, ios, npy, naky, my, jj, kk, jjj&
         , kkk, ierr
      real(real64) :: e(max_orbitals), gama(max_orbitals)
      real(real64) :: con, fkk, ey, pzy, dnorm
      real(real64), pointer :: pa(:), qa(:), ra(:)
      character  :: nh(max_orbitals)*2
      character :: g92rwf*6
!-----------------------------------------------
!
!  Common relevant for the initial state
!
!
      ntestl = 0
      ntest = max0(ntestl,ntestg)
      ntest = 0

!
!   Write entry message
!
      write (6, *) 'Loading Radial WaveFunction File for initial state...'
!
!   Open the radial wave function file
!
      j = index(name,' ')
      open(unit=21, file=name(1:j-1)//'.w', form='UNFORMATTED', status='OLD', &
         position='asis')
!
!   Save NAK, NP and NH
!
      nak(:nwii) = nakii(:nwii)
      np(:nwii) = npii(:nwii)
      nh(:nwii) = nhii(:nwii)
!
!   Allocate storage to orbital arrays
!
      call alloc (pfii, max_radial_points,nwii, 'PFII', 'LODRWII')
      call alloc (qfii, max_radial_points,nwii, 'QFII', 'LODRWII')

!
      con = nuclear_charge/speed_of_light_au
      con = con*con
!
      do j = 1, nwii
         pfii(:max_radial_points,j) = 0.0_real64
         qfii(:max_radial_points,j) = 0.0_real64
!
         k = abs(nak(j))
         if (nuclear_model_parameter_count /= 0) cycle
         fkk = dble(k*k)
         if (fkk >= con) then
            gama(j) = sqrt(fkk - con)
         else
            write (6, *) 'LODRWF: Imaginary gamma parameter'
            write (6, *) ' for ', np(j), nh(j), ' orbital; the'
            write (6, *) ' point model for the nucleus'
            write (6, *) ' is inappropriate for Z > ', speed_of_light_au, '.'
            stop
         endif
!
      end do
!
!   Read orbital information from Read Orbitals File;
!
      nwin = 0
      read (21, iostat=ios) g92rwf
      if (ios/=0 .or. g92rwf/='G92RWF') then
         write (6, *) 'This is not a Radial WaveFunction File;'
         close(21)
      endif

      if (ntest >= 100) then
         write (*, *) '******************'
         write (*, *) ' Entering lodrwfi'
         write (*, *) '******************'
      endif

    3 continue
      read (21, iostat=ios) npy, naky, ey, my
      if (ios == 0) then
         call alloc (pa, my, 'PA', 'LODRWFII')
         call alloc (qa, my, 'QA', 'LODRWFII')
         call alloc (ra, my, 'RA', 'LODRWFII')
         read (21) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
         read (21) (ra(i),i=1,my)
!
!    Orbital order as defined in kapdata
!
         jj = 0
         do k = 1, kamax
            if (k > 1) jj = nshlii(k-1) + jj
            do j = 1, nshlii(k)
               kk = nshlpii(k,j)
               if (npy/=np(kk) .or. naky/=nak(kk)) cycle
               jjj = jj + j
               pzii(jjj) = pzy
               eii(jjj) = ey
               nakii(jjj) = nak(kk)
               npii(jjj) = np(kk)
               nhii(jjj) = nh(kk)
               gamaii(jjj) = gama(kk)
               call intrpqi (pa, qa, my, ra, jjj, dnorm)
               if (ntest >= 100) write (*, 301) npii(jjj), nhii(jjj), eii(jjj)&
                  , dnorm
               if (ntest > 1000) then
                  write (*, *) 'PF              QF             RA'
                  do kkk = 1, mfii(jjj)
                     write (*, *) pfii(kkk,jjj), qfii(kkk,jjj), ra(kkk)
                  end do
               endif
               nwin = nwin + 1
            end do
         end do
         call dalloc (pa, 'PA', 'LODRWII')
         call dalloc (qa, 'QA', 'LODRWII')
         call dalloc (ra, 'RA', 'LODRWII')
         go to 3
      endif
      if (ldbpr(3)) write (99, *) ' orbitals renormalised;'
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nwii) then
         write (6, *) 'LODRWF: All required orbitals not'
         write (6, *) ' found.'
         ierr = 1
         go to 5
      endif
!
      write (6, *) ' ... load complete;'
!
    5 continue
      close(21)
      if (ntest >= 100) then
         write (*, *) 'Sorted order should be the same as from kapdat'
         do j = 1, nwii
            write (*, 301) npii(j), nhii(j), eii(j), dnorm
         end do
         write (*, *)
         write (*, *) '*****************'
         write (*, *) ' Leaving lodrwfi'
         write (*, *) '*****************'
      endif

      return
!
  301 format(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)
      return
!
      end subroutine lodrwfi
