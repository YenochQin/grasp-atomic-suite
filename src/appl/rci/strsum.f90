!***********************************************************************
!                                                                      *
      subroutine strsum
!                                                                      *
!   Generates the first part of  grasp92.sum  (on stream 24).          *
!                                                                      *
!   Call(s) to: [LIB92] CALEN, convrt2.                                *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 09 Dec 1992   *
!   Modified by Xinghong He               Last revision: 22 Dec 1997   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use decide_c
      use def_c
      use grid_c
      use npar_c
      use npot_c
      use orb_c
      use prnt_c
      use wave_c
      use wfac_c
      use iounit_c
      use hblock_c
      use iccu_c,          only: iccutblk
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt2_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: lenth, nb, iccut, i, iend, ibeg
      character :: record*256, cdata*26, clevel*2
!-----------------------------------------------
!     POINTER (pncfblk, ncfblk(0:*))
!
!     POINTER (piccutblk, iccutblk(1))

!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!   Write out the basic dimensions of the electron cloud
!
      write (24, *)
      call convrt2 (electron_count, record, lenth, 'strsum.nelec')
      write (24, *) 'There are '//record(1:lenth)//' electrons in the cloud'
      call convrt2 (csf_count, record, lenth, 'strsum.ncf')
      write (24, *) ' in '//record(1:lenth)//' relativistic CSFs'
      call convrt2 (orbital_count, record, lenth, 'strsum.nw')
      write (24, *) ' based on '//record(1:lenth)//' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      if (lfordr) then
         write (24, *)
         do nb = 1, symmetry_block_count
            iccut = iccutblk(nb)
            call convrt2 (iccut, record, lenth, 'strsum.icccut')
            write (24, *) ' CSFs 1--'//record(1:lenth)//' constitute'//&
               ' the zero-order space;  nb = ', nb, ' ncf = ', csf_counts_by_block(nb)
         end do
      endif
!
!   Write out the nuclear parameters
!
      write (24, *)
      write (24, 300) nuclear_charge
      if (emn == 0.0d00) then
         write (24, *) ' the nucleus is stationary;'
      else
         write (24, 301) emn
      endif
      if (nuclear_model_parameter_count == 2) then
         write (24, *) ' Fermi nucleus:'
         write (24, 302) nuclear_model_parameters(1), nuclear_model_parameters(2)
         call convrt2 (first_grid_index_beyond_fermi_c, record, lenth, 'strsum.nnuc')
         write (24, *) ' there are '//record(1:lenth)//&
            ' tabulation points in the nucleus.'
      else
         write (24, *) ' point nucleus.'
      endif
!
!   Write out the physical effects specifications
!
      write (24, *)
      write (24, 303) speed_of_light_au
!
      write (24, *)
      if (ltrans .or. lvp .or. lnms .or. lsms) then
         write (24, *) 'To H (Dirac Coulomb) is added'
         if (ltrans) write (24, 304) wfact
         if (lvp) write (24, *) ' H (Vacuum Polarisation);'
         if (lnms) write (24, *) ' H (Normal Mass Shift);'
         if (lsms) write (24, *) ' H (Specific Mass Shift);'
         write (24, *) ' the total will be diagonalised.'
      else
         write (24, *) 'H (Dirac Coulomb) will be diagonalised by itself.'
      endif
!
      if (lse) then
         write (24, *) &
            'Diagonal contributions from H (Self Energy) will be estimated'
         write (24, *) ' from a screened hydrogenic approximation.'
      endif
!
!   Write out the parameters of the radial grid
!
      write (24, *)
      if (asymptotic_grid_step == 0.0d00) then
         write (24, 305) radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (24, 306) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      endif
      write (24, 307) radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
!
!   Write out the orbital properties
!
      write (24, *)
      write (24, *) 'Subshell radial wavefunction summary:'
      write (24, *)
      write (24, 308)
      write (24, *)
      do i = 1, orbital_count
         write (24, 309) orbital_principal_n(i), orbital_symmetry_labels(i), orbital_energy_parameters(i), &
            & orbital_origin_parameters(i), orbital_origin_exponents(i), radial_large_components(2,i), &
            & radial_small_components(2,i), &
            orbital_last_grid_index(i)
      end do
!
!   Write the list of eigenpair indices
!
      write (24, *)
      if (nvec == 1) then
         call convrt2 (ivec(1), record, lenth, 'strsum.ivec(1)')
         write (24, *) 'Level '//record(1:lenth)//' will be computed.'
      else
         call convrt2 (nvec, record, lenth, 'strsum.nvec')
         write (24, *) record(1:lenth)//' levels will be computed;'
         record(1:20) = ' their indices are: '
         iend = 20
         do i = 1, nvec
            ibeg = iend + 1
            call convrt2 (ivec(i), clevel, lenth, 'strsum.ivec(i)')
            if (i /= nvec) then
               iend = ibeg + lenth + 1
               record(ibeg:iend) = clevel(1:lenth)//', '
            else
               iend = ibeg + lenth
               record(ibeg:iend) = clevel(1:lenth)//'.'
            endif
            if (iend < 120) cycle
            write (24, *) record(1:iend)
            record(1:2) = '  '
            iend = 2
         end do
         if (iend /= 2) write (24, *) record(1:iend)
      endif
!
      return
!
  300 format('The atomic number is ',1f14.10,';')
  301 format(' the mass of the nucleus is ',1p,d19.12,' electron masses;')
  302 format('  c =',1p,1d19.12,' Bohr radii,'/,'  a =',1d19.12,' Bohr radii;')
  303 format('Speed of light = ',1p,d19.12,' atomic units.')
  304 format('  H (Transverse) --- factor multiplying the',&
         ' photon frequency: ',1p,d15.8,';')
  305 format('Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',' I = 1, ..., N;'/,/,&
         ' RNT  = ',1p,d19.12,' Bohr radii;'/,' H    = ',d19.12,' Bohr radii;'/&
         ,' N    = ',1i4,';')
  306 format('Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',&
         ' I = 1, ..., N;'/,/,' RNT  = ',1p,d19.12,' Bohr radii;'/,' H    = ',d&
         19.12,' Bohr radii;'/,' HP   = ',d19.12,' Bohr radii;'/,' N    = ',1i4&
         ,';')
  307 format(' R(1) = ',1p,1d19.12,' Bohr radii;'/,' R(2) = ',1d19.12,&
         ' Bohr radii;'/,' R(N) = ',1d19.12,' Bohr radii.')
  308 format('Subshell',6x,'e',13x,'p0',5x,'gamma',5x,'P(2)',7x,'Q(2)',4x,'MTP'&
         )
  309 format(1x,i2,a2,1x,1p,d17.10,1p,d11.3,0p,f6.2,1p,2(d11.3),i5)
      return
!
      end subroutine strsum
