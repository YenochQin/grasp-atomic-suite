!***********************************************************************
!                                                                      *
      subroutine lodres
!                                                                      *
!   Loads the data from the  .res  file. A number of checks are made   *
!   to ensure correctness and consistency.                             *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, GETYN, SETQIC.                 *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 15 Oct 1992   *
!   Block version by Xinghong He          Last revision:  1 Jun 1998   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use parameter_def,   only: nnnp
      use memory_man
      use decide_c
      use def_c, only: nelecr, nelec, z, emn, c, accy
      use grid_c
      use npar_c
      use npot_c, only: zz, nnuc
      use orb_c
      use wave_c
      use wfac_c
      use where_c
      use hblock_c
      use iccu_c
      use iounit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use getyn_i
      use setqic_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ncfres, nwres, nblockres, i, np10, j
      logical :: yes
!-----------------------------------------------
!
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))
!
!-----------------------------------------------------------------------
      write (6, *) 'Calling lodres ...'
!
!   Read the basic parameters of the electron cloud; check these
!   against those deduced from the  .csl file
!
      read (imcdf) nelecr, ncfres, nwres, nblockres

      if (nelecr/=nelec .or. ncfres/=ncf .or. nwres/=nw .or. nblockres/=nblock&
         ) stop 'lodres: NELEC/NCF/NW does not match'
!
!   Read the nuclear parameters
!
      read (imcdf) z, emn
      read (imcdf) nparm, (parm(i),i=1,nparm)
      read (imcdf) n, (zz(i),i=1,n), nnuc

      if (n > nnnp) stop 'lodres: N greater than NNNP'
!
!   Read the physical effects specifications
!   iccutblk() is now an array of length nblock.
!
      read (imcdf) c, lfordr, (iccutblk(i),i=1,nblock), ltrans, wfact, lvp, &
         lnms, lsms
!
!   Read the remaining parameters controlling the radial grid and the
!   grid arrays
!
      np10 = n + 10
      read (imcdf) rnt, h, hp, (r(i),i=1,np10), (rp(i),i=1,np10), (rpor(i),i=1,&
         np10)
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = h**6
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Allocate storage for the radial wavefunction arrays
!
      call alloc (pf, nnnp,nw, 'PF', 'LODMIX')
      call alloc (qf, nnnp,nw, 'QF', 'LODMIX')
!
!   Read the orbital wavefunctions and the associated arrays
!
      do j = 1, nw
         read (imcdf) e(j), gama(j), pz(j), mf(j)
         read (imcdf) (pf(i,j),i=1,mf(j)), (qf(i,j),i=1,mf(j))
      end do
!
!   Determine if the self-energy contribution is to be estimated
!
      write (istde, *) 'Estimate contributions from the self-energy?'
      lse = getyn()
      if (lse) then
        write(734,'(a)') 'y            ! Estimate contributions from the self-energy?'
      else
        write(734,'(a)') 'n            ! Estimate contributions from the self-energy?'
      end if

      return
      end subroutine lodres
