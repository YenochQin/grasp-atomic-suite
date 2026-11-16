!***********************************************************************
!                                                                      *
      subroutine lodrwf(iierr)
!                                                                      *
!   This subroutine loads  radial wavefunctions from the  .rwf  file   *
!   and performs some related setup.                                   *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, DALLOC, INTRPQ, ORTHSC.                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
!   Block version by Xinghong He          Last revision: 27 May 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:24:43   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use memory_man
      use debug_c
      use def_c,           only: c, z
      use grid_c
      use npar_c
      use orb_c
      use wave_c,          only: pz, pf,qf
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use intrpq_i
      use orthsc_i
      implicit none
!-----------------------------------------------
!   G l o b a l   P a r a m e t e r s
!-----------------------------------------------
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out) :: iierr
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: j, i, k, nwin, ios, npy, naky, my
      integer :: ierr
      real(double) :: con, fkk, ey, pzy, dnorm
      real(double), dimension(:), pointer :: pa, qa, ra
!-----------------------------------------------
!
!
!   Write entry message
!
      write (6, *) 'Loading Radial WaveFunction File ...'
!
!   Allocate storage to orbital arrays
!
      call alloc(pf, nnnp, nw, 'PF', 'LODRWF')
      call alloc(qf, nnnp, nw, 'QF', 'LODRWF')

!
!   Setup: (1) Orbital arrays to zero
!          (2) Array E to -1 (no orbitals estimated)
!          (3) Parameters GAMMA for each orbital
!
      con = z/c
      con = con*con
!
      do j = 1, nw
!
         pf(:n,j) = 0.0d00
         qf(:n,j) = 0.0d00
!
         e(j) = -1.0d00
!
         k = abs(nak(j))
         if (nparm > 0) then
            gama(j) = dble(k)
         else if (nparm == 0) then
            fkk = dble(k*k)
            if (fkk >= con) then
               gama(j) = sqrt(fkk - con)
            else
               !WRITE (istde,*) 'LODRWF: Imaginary gamma parameter'
               !WRITE (istde,*) ' for ',NP(J),NH(J),' orbital; the'
               !WRITE (istde,*) ' point model for the nucleus'
               !WRITE (istde,*) ' is inappropriate for Z > ',C,'.'
               stop 'lodrwf: Inappropriate gamma'
            endif
         endif
!
      end do
!
!   Read orbital information from Read Orbitals File; write summary
!   to  .dbg  file if option set
!
      if (ldbpr(3)) write (99, 300)
      nwin = 0
    3 continue
      read (23, iostat=ios) npy, naky, ey, my

      if (ios == 0) then
         call alloc (pa,my, 'PA', 'LODRWF')
         call alloc (qa,my, 'QA', 'LODRWF')
         call alloc (ra,my, 'RA', 'LODRWF')

         read (23) pzy, (pa(i),i=1,my), (qa(i),i=1,my)
         read (23) (ra(i),i=1,my)

         do j = 1, nw
            if (.not.(e(j)<0.0d00 .and. npy==np(j) .and. naky==nak(j))) cycle
            pz(j) = pzy
            e(j) = ey
            call intrpq(pa, qa, my, ra, j, dnorm)
            if (ldbpr(3)) write (99, 301) np(j), nh(j), e(j), dnorm
            nwin = nwin + 1
         end do
         call dalloc (pa, 'PA', 'LODRWF' )
         call dalloc (qa, 'QA', 'LODRWF')
         call dalloc (ra, 'RA', 'LODRWF')
         go to 3
      endif
      if (ldbpr(3)) write (99, *) ' orbitals renormalised;'
!
!   Stop with an error message if all orbitals are not known
!
      if (nwin < nw) then
         iierr = 1
         go to 5
      endif
!
!   Schmidt orthogonalise the orbitals
!
      call orthsc
      if (ldbpr(3)) write (99, *) ' orbitals orthogonalised and renormalised;'
!
      iierr = 0
    5 continue
      return
!
  300 format(/,'From SUBROUTINE LODRWF:'/,' Orbital',8x,'Eigenvalue',19x,'Norm')
  301 format(2x,i2,a2,4x,1p,1d22.15,4x,1d22.15)
      return
!
      end subroutine lodrwf
