!***********************************************************************
!                                                                      *
      subroutine strsum(name, inpci, ilbl)
      use calculation_file_units, only: lsj_transition_output_unit, summary_output_unit
!                                                                      *
!   Generates the first part of  oscl92.sum  (on stream 24).           *
!                                                                      *
!   Call(s) to: [LIB92]: CALEN, CONVRT, WGHTD5.                        *
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 28 Dec 1992   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use biorthogonal_orbital_sets, only: max_orbitals, nhff, nhii, eff, gamaff, eii, gamaii, ncfff, nwff, ncfii, nwii, npff, &
          nakff, npii, nakii, nklff, nkjff, nklii, nkjii, nnff, nnii
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use first_order_csf_cutoffs, only: iccut
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use orbital_csf_basis, only: max_orbitals
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use level_symmetry_metadata, only: iaspar, iasparff, iasparii, iatjpo, iatjpoff, iatjpoii, jpgg, nblk0
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use transition_lsj_labels, only: rlev_ener_1, rlev_ener_2, string_csf1, string_csf2, nvectoti, nvectotf, &
          iopen_status1, iopen_status2, lev_pos_1, lev_j_1, lev_par_1, lev_pos_2, lev_j_2, lev_par_2
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: inpci
      integer :: ilbl
      character, intent(in) :: name(2)*24
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i, j
      character(len=:), allocatable :: lsj_transition_filename
      character :: record*15, time_text*8, date_text*8
!-----------------------------------------------


      i = index(name(1)," ")
      j = index(name(2)," ")
      if(ilbl == 0) then
         if (inpci == 0) then
            open(newunit=summary_output_unit,file=name(1)(1:i-1)//"."//name(2)(1:j-1)//".ct",form=&
            "FORMATTED", status="UNKNOWN",position="asis")
         else
            open(newunit=summary_output_unit,file=name(1)(1:i-1)//"."//name(2)(1:j-1)//".t",form=&
            "FORMATTED",status="UNKNOWN",position="asis")
         end if
      else if(ilbl == 1) then
         if(iopen_status1==0 .and. iopen_status2 ==0) then
            if (inpci == 0) then
               lsj_transition_filename = name(1)(1:i-1)//"."//name(2)(1:j-1)//".ct.lsj"
            else
               lsj_transition_filename = name(1)(1:i-1)//"."//name(2)(1:j-1)//".t.lsj"
            end if
            if (lsj_transition_output_unit == -1) then
               open(newunit=lsj_transition_output_unit, file=lsj_transition_filename, &
                  form="FORMATTED", status="UNKNOWN")
            else
               open(unit=lsj_transition_output_unit, file=lsj_transition_filename, &
                  form="FORMATTED", status="UNKNOWN")
            end if
         end if
      end if
!
!   Get the date and time of day; make this information the
!   header of the summary file
!
!ww      CALL CALEN (CTIME,CDATE)
!ww      WRITE (24,*) 'OSCL92 run at ',CTIME,' on ',CDATE,'.'
!
!   Write out the basic dimensions of the initial state electron cloud
!
!ww      WRITE (24,*)
!ww      CALL CONVRT (NELECII,RECORD,LENTH)
!ww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
!ww     :           //' electrons in the initial state cloud'
!ww      CALL CONVRT (NCFII,RECORD,LENTH)
!ww      WRITE (24,*) ' in '//RECORD(1:LENTH)
!ww     :           //' relativistic CSFs'
!ww      CALL CONVRT (NWII,RECORD,LENTH)
!ww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
!ww     :           //' relativistic subshells.'
!
!   Write out the basic dimensions of the final state electron cloud
!
!ww      WRITE (24,*)
!ww      CALL CONVRT (NELECFF,RECORD,LENTH)
!ww      WRITE (24,*) 'There are '//RECORD(1:LENTH)
!ww     :           //' electrons in the final state cloud'
!ww      CALL CONVRT (NCFFF,RECORD,LENTH)
!ww      WRITE (24,*) ' in '//RECORD(1:LENTH)
!ww     :           //' relativistic CSFs'
!ww      CALL CONVRT (NWFF,RECORD,LENTH)
!ww      WRITE (24,*) ' based on '//RECORD(1:LENTH)
!ww     :           //' relativistic subshells.'
!
!   If the CSFs are not treated uniformly, write out an
!   informative message
!
!ww      IF (LFORDR) THEN
!ww         WRITE (24,*)
!ww         CALL CONVRT (ICCUT,RECORD,LENTH)
!ww         WRITE (24,*) ' CSFs 1--'//RECORD(1:LENTH)//' constitute'
!ww     :              //' the zero-order space;'
!ww      ENDIF
!
!   Write out the nuclear parameters
!
!ww      WRITE (24,*)
!ww      WRITE (24,300) Z
!ww      IF (NPARM .EQ. 2) THEN
!ww         WRITE (24,*) 'Fermi nucleus:'
!ww         WRITE (24,301) PARM(1),PARM(2)
!ww      ELSE
!ww         WRITE (24,*) ' point nucleus.'
!ww      ENDIF
!
!   Write out the physical effects specifications
!
!ww      WRITE (24,*)
!ww      WRITE (24,305) C
!
!   Write out the parameters of the radial grid
!
!ww      WRITE (24,*)
!ww      IF (HP .EQ. 0.0D 00) THEN
!ww         WRITE (24,306) RNT,H,N
!ww      ELSE
!ww         WRITE (24,307) RNT,H,HP,N
!ww      ENDIF
!ww      WRITE (24,308) R(1),R(2),R(N)
!
!   Write out the orbital properties
!
!ww      WRITE (24,*)
!ww      WRITE (24,*) 'Initial state subshell radial wavefunction summary:'
!ww      WRITE (24,*)
!ww      WRITE (24,309)
!ww      WRITE (24,*)
!ww      DO 1 I = 1,NWII
!ww         WRITE (24,310) NPII(I),NHII(I),EII(I),PZII(I),
!ww     :                  GAMAII(I),PFII(2,I),QFII(2,I),MFII(I)
!ww    1 CONTINUE
!
!ww      WRITE (24,*)
!ww      WRITE (24,*) 'Final state subshell radial wavefunction summary:'
!ww      WRITE (24,*)
!ww      WRITE (24,309)
!ww      WRITE (24,*)
!ww      DO 2 I = 1,NWFF
!ww         WRITE (24,310) NPFF(I),NHFF(I),EFF(I),PZFF(I),
!ww     :                  GAMAFF(I),PFFF(2,I),QFFF(2,I),MFFF(I)
!ww    2 CONTINUE
!
!   Write the list of eigenpair indices for the initial state
!
!     WRITE (24,*)
!     CALL ENGOUT1 (EAVII,EVALII,IATJPOII,IASPARII,IVECII,NVECII,3,1)
!
!   Write the list of eigenpair indices for the final state
!
!     WRITE (24,*)
!     CALL ENGOUT1 (EAVFF,EVALFF,IATJPOFF,IASPARFF,IVECFF,NVECFF,3,2)
!
      return
!

      return
!
      end subroutine strsum
