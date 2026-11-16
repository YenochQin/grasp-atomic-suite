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
!...Translated by Gediminas Gaigalas 11/18/19
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
      use syma_c
      use npar_c
      use eigv_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use convrt2_i
      use engout_i
      implicit none
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer   :: lenth, nb, iccut, i, iend, ibeg
      character :: record*256, cdata*26, clevel*2
!-----------------------------------------------
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
         call convrt2 (iccut, record, lenth, 'strsum.icccut')
         write (24, *) ' CSFs 1--'//record(1:lenth)//' constitute'//   &
              ' the zero-order space;'
      endif
!
!   Write out the nuclear parameters
!
      write (24, *)
      write (24, 300) nuclear_charge
      if (nuclear_model_parameter_count == 2) then
         write (24, *) ' Fermi nucleus:'
         write (24, 301) nuclear_model_parameters(1), nuclear_model_parameters(2)
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
      write (24, 305) speed_of_light_au
!
!   Write out the parameters of the radial grid
!
      write (24, *)
      if (asymptotic_grid_step == 0.0d00) then
         write (24, 306) radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (24, 307) radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      endif
      write (24, 308) radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
!
!   Write out the orbital properties
!
      write (24, *)
      write (24, *) 'Subshell radial wavefunction summary:'
      write (24, *)
      write (24, 309)
      write (24, *)
      do i = 1, orbital_count
         write (24, 310) orbital_principal_n(i), orbital_symmetry_labels(i), orbital_energy_parameters(i), &
            & orbital_origin_parameters(i), orbital_origin_exponents(i), radial_large_components(2,i), &
            & radial_small_components(2,i), &
            orbital_last_grid_index(i)
      end do
!
!   Write the list of eigenpair indices
!
      write (24,*)
      call engout (hamiltonian_energy_shift,shifted_eigenvalues,iatjpo,iaspar,ivec,nvec,3)
!
      return
!
  300 format ('The atomic number is ',1f14.10,';')
  301 format ('  c =',1p,1d19.12,' Bohr radii,'                        &
             /'  a =',   1d19.12,' Bohr radii;')
  305 format ('Speed of light = ',1pd19.12,' atomic units.')
  306 format ( 'Radial grid: R(I) = RNT*(exp((I-1)*H)-1),',            &
               ' I = 1, ..., N;'                                       &
             //' RNT  = ',1p,d19.12,' Bohr radii;'                     &
              /' H    = ',   d19.12,' Bohr radii;'                     &
              /' N    = ',1i4,';')
  307 format ( 'Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,',   &
               ' I = 1, ..., N;'                                       &
             //' RNT  = ',1p,d19.12,' Bohr radii;'                     &
              /' H    = ',   d19.12,' Bohr radii;'                     &
              /' HP   = ',   d19.12,' Bohr radii;'                     &
              /' N    = ',1i4,';')
  308 format ( ' R(1) = ',1p,1d19.12,' Bohr radii;'                    &
              /' R(2) = ',   1d19.12,' Bohr radii;'                    &
              /' R(N) = ',   1d19.12,' Bohr radii.')
  309 format (' Subshell',11x,'e',20x,'p0',18x,                        &
              'gamma',19x,'P(2)',18x,'Q(2)',10x,'MTP')
  310 format (3x,1i2,1a2,1x,1p,5(3x,1d19.12),3x,1i3)
!
      end subroutine strsum
