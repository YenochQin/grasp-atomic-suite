!***********************************************************************
      subroutine strsum

!   Generates the first part of  rscf92.sum  (on stream 24).
!
!   Call(s) to: [LIB92] CALEN, CONVRT.
!
!   Written by Farid A. Parpia            Last revision: 26 Sep 1993
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: electron_count, emn, iccmin, ioncty, max_radial_points, ncmin, nuclear_charge, &
          speed_of_light_au, weight
      use first_order_csf_cutoffs, only: iccut
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use rmcdhf_mcp_options, only: diag, lfordr
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_potential_state, only: first_grid_index_beyond_fermi_c, max_radial_points, negative_radius_times_potential
      use orbital_csf_basis, only: csf_count, max_orbitals, orbital_count
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use breit_frequency_scale, only: wfact
!      USE iccu_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use integer_format, only: convrt
      implicit none (type, external)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: text_length, iend, i, ibeg
      character :: record*256, cdata*26, time_text*8, date_text*8, number_text*4
!-----------------------------------------------
!
!     POINTER (PWEIGH,WEIGHT(1))
!     POINTER (PCCMIN,ICCMIN(1))
!     POINTER (PNTRPF,PF(NNNP,1))
!     POINTER (PNTRQF,QF(NNNP,1))
!
!   Both the nuclear charge and the number of electrons are
!   known at this point; load IONCTY with the ionicity
!
      ioncty = nint(nuclear_charge) - electron_count
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!      CALL CALEN (CTIME, CDATE)
!      WRITE (24,*) 'RSCF92 run at ',CTIME,' on ',CDATE,'.'
!
!   Write out the basic dimensions of the electron cloud
!
      write (24, *)
      call convrt (electron_count, record, text_length)
      write (24, *) "There are "//record(1:text_length)//" electrons in the cloud"
      call convrt (csf_count, record, text_length)
      write (24, *) " in "//record(1:text_length)//" relativistic CSFs"
      call convrt (orbital_count, record, text_length)
      write (24, *) " based on "//record(1:text_length)//" relativistic subshells."
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
      if (lfordr) then
         write (24, *)
         call convrt (iccut, record, text_length)
         write (24, *) " CSFs 1--"//record(1:text_length)//" constitute"//&
            " the zero-order space;"
      end if
!
!   Write out the nuclear parameters
!
      write (24, *)
      write (24, '("The atomic number is ",1f14.10,";")') nuclear_charge

      if (emn == 0.0_real64) then
         write (24, *) " the nucleus is stationary;"
      else
         write (24, '(" the mass of the nucleus is ",1p,d19.12," electron masses;")') emn
      end if
      if (nuclear_model_parameter_count == 2) then
         write (24, *) " Fermi nucleus:"
         write (24, '("  c =",1p,1d19.12," Bohr radii,"/,"  a =",1d19.12," Bohr radii;")') &
            nuclear_model_parameters(1), nuclear_model_parameters(2)
         call convrt (first_grid_index_beyond_fermi_c, record, text_length)
         write (24, *) " there are "//record(1:text_length)//&
            " tabulation points in the nucleus."
      else
         write (24, *) " point nucleus."
      end if
!
!   Write out the physical effects specifications
!
      write (24, *)
      write (24, '("Speed of light = ",3p,d19.12," atomic units.")') speed_of_light_au
!
!   Write out the parameters of the radial grid
!
      write (24, *)
      if (asymptotic_grid_step == 0.0_real64) then
         write (24, '("Radial grid: R(I) = RNT*(exp((I-1)*H)-1),"," I = 1, ..., N;"/,/," RNT  = ",1p,d19.12,&
            &" Bohr radii;"/," H    = ",d19.12," Bohr radii;"/," N    = ",1i4,";")') &
            radial_grid_scale, transformed_grid_step, radial_point_count
      else
         write (24, '("Radial grid: ln(R(I)/RNT+1)+(H/HP)*R(I) = (I-1)*H,"," I = 1, ..., N;"/,/," RNT  = ",1p,&
            &d19.12," Bohr radii;"/," H    = ",d19.12," Bohr radii;"/," HP   = ",d19.12," Bohr radii;"/," N    = ",&
            &1i4,";")') radial_grid_scale, transformed_grid_step, asymptotic_grid_step, radial_point_count
      end if
      write (24, '(" R(1) = ",1p,1d19.12," Bohr radii;"/," R(2) = ",1d19.12," Bohr radii;"/," R(N) = ",1d19.12,&
         &" Bohr radii.")') radial_grid(1), radial_grid(2), radial_grid(radial_point_count)
      write (24, *)
!
!  (E)AL calculation, returns here
!
      if (ncmin == 0) then
         write (24, *) "(E)AL calculation."
         return
      end if
!
!  Info exclusively for EOL calculations
!
      if (ncmin == 1) then
         write (24, *) "OL calculation."
         call convrt (iccmin(1), record, text_length)
         write (24, *) "Level "//record(1:text_length)//" will be optimised."
      else
         write (24, *) "EOL calculation."
         call convrt (ncmin, record, text_length)
         write (24, *) record(1:text_length)//" levels will be optimised;"
         record(1:20) = " their indices are: "
         iend = 20
         do i = 1, ncmin
            ibeg = iend + 1
            call convrt (iccmin(i), number_text, text_length)
            if (i /= ncmin) then
               iend = ibeg + text_length + 1
               record(ibeg:iend) = number_text(1:text_length)//", "
            else
               iend = ibeg + text_length
               record(ibeg:iend) = number_text(1:text_length)//"."
            end if
            if (iend < 120) cycle
            write (24, *) record(1:iend)
            record(1:2) = "  "
            iend = 2
         end do
         if (iend /= 2) write (24, *) record(1:iend)
         if (weight(1) == (-1.0_real64)) then
            write (24, *) "Each is assigned its statistical weight;"
         else if (weight(1) == (-2.0_real64)) then
            write (24, *) "All levels are weighted equally;"
         else
            write (24, *) " weighted as follows:"
            write (24, *) (weight(i),i=1,ncmin)
         end if
      end if


      return
      end subroutine strsum
