!***********************************************************************
!                                                                      *
      subroutine setham (myid, nprocs, jblock, elsto,icstrt, nelmntt, &
                         atwinv,slf_en)
!                                                                      *
!   Sets up the Hamiltonian matrix and determines the average energy.  *
!
!   Serial I/O moved out; able to run on single/multiple processors
!   For this purpose a new common /setham_to_genmat2/ is created
!                                                                      *
!   Call(s) to: [LIB92]: ALCBUF, CONVRT, DALLOC, ICHOP, RKCO, TNSRJJ.  *
!               [RCI92]: BRINT, IABINT, KEINT, RKINTC, VINT, VPINT.    *
!                                                                      *
!   Written by Farid A Parpia             Last revision: 30 Oct 1992   *
!   Block version by Xinghong He          Last revision: 15 Jun 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use grasp_parameters,   only: max_radial_points, orbital_index_radix
      use array_memory_management, only: alloc, dalloc
!cjb iso_fortran_env :: output_unit, error_unit
      use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
!-----------------------------------------------
!   C O M M O N    B l o c k s
!-----------------------------------------------
      use breit_core_occupations, only: icore
      use breit_integral_tables, only: ndtpa, ntpi, indtp1, indtp2, indtp3, indtp4, indtp5, indtp6, &
          valtp1, valtp2, valtp3, valtp4, valtp5, valtp6, first
      use angular_coefficient_buffer_state, only: angular_buffer_capacity, angular_coefficient_count, &
          angular_coefficient_labels, angular_coefficients
      use coulomb_integral_table, only: frstco, indoei, ncoei, ndcoea, valoei
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use debug_output_settings, only: cutoff, ibug1, ibug2, ibug3, ibug4, ibug5, ibug6, ldbpa, ldbpg, ldbpr
      use atomic_calculation_state, only: emn, max_radial_points
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use csf_block_cutoffs, only: iccut, iccutblk
      use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
          radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
      use sparse_hamiltonian_storage, only: stored_matrix_element_count, hamiltonian_values, stored_segment_end_offsets, &
          matrix_partner_indices, ntpitmp, nelmnttmp, ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp, ncoretmp, &
          nvpitmp, nkeitmp, nvintitmp, ncftmp, cutofftmp
      use kinetic_energy_integral_table, only: frstki, indkei, ndkea, nkei, valkei
      use nuclear_charge_distribution, only: max_radial_points, zdist
      use orbital_csf_basis, only: csf_count,  orbital_count
      use selected_level_indices, only: ivec, ivecff, ivecii, nvec, nvecff, nvecii, nvecmx, nvecmxff, nvecmxii, pnivecii
      use csf_coupling_data, only: csf_coupling_quantum_codes, csf_subshell_quantum_codes, jcupar, jqsar
      use breit_bessel_cache_keys, only: keep
      use radial_quadrature_workspace, only: mtp, radial_workspace_capacity, ta, tb
      use vinti_integral_table, only: frstvi, indtei, ndvin, nvinti, valtei
      use vacuum_polarization_integral_table, only: frstvp, indvpi, ndvpa, nvpi, valvpi
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
      use distributed_slater_integrals, only: indteirk, valteirk, rk_local_count
      use rci_prerun_settings, only: coeffcut1, coeffcut2, iprerun, ncsfpre
      use where_c, only: imcdf
      use scf_previous_eigenvectors, only: evec1
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angular_coefficient_buffer, only: alcbuf
      use iabint_i, only: iabint
      use brint1_i, only: brint1
      use brint2_i, only: brint2
      use brint3_i, only: brint3
      use brint4_i, only: brint4
      use brint5_i, only: brint5
      use brint6_i, only: brint6
      use integer_format, only: convrt
      use subshell_occupation_class, only: ichop
      use keint_i, only: keint
      use rkintc_i, only: rkintc
      use vint_i, only: vint
      use vpint_i, only: vpint
      use one_electron_scalar_coefficients, only: onescalar
      use two_electron_angular_coefficients, only: rkco_gg
      use closed_shell_mcp, only: cord
      implicit none (type, external)
      external :: breid
!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer(int64)       :: nelmntt
      integer             :: jblock, icstrt
      integer, intent(in) :: myid, nprocs
      real(real64)        :: elsto, atwinv
      real(real64) :: slf_en(*)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      real(real64) :: tshell(max_orbitals)
      real(real64) :: tgrl1, tgrl2, tegral

      integer, parameter :: key = orbital_index_radix
!
!     Matrix elements smaller than CUTOFF are not accumulated
!
!cjb cutoff is use associated and cannot be redeclared
!     REAL(DOUBLE), PARAMETER :: CUTOFF = 1.0D-20
!cjb CUTOFF = 1.0D-12 below
!cjb
      integer :: ipi, ipj, inc1, inc2, kt, ipt, incor, ncoec, nctec, &
                 i, j, nmcbp, ncore, ic, nelc, irstart, ir, ia, ib,  &
                 itype, iia
       real(real64) :: elemnt, precoeff, tcoeff, vcoeff, contr
!-----------------------------------------------------------------------
      stored_matrix_element_count = nelmntt
!cjb
!cjb  CUTOFF = 1.0D-12
      cutoff = 1.0e-12_real64

      atwinv = 1.0_real64/emn

      if (iprerun == 2) then
         do ipi = 1,nvec
            do ipj = 1,csf_count
               write (*,*) ipi,ipj,evec1(ipj+(ipi-1)*csf_count)
            end do
         end do
      end if
!
!   Allocate storage to arrays in COMMON/BUFFER/; these are
!   used for the Coulomb and transverse two-electron integrals
!
      call alcbuf (1)

!     ...Locals
      call alloc (hamiltonian_values, csf_count, "EMT","SETHAM" )
      call alloc (matrix_partner_indices, csf_count, "IROW", "SETHAM")
!
      inc1 = 1
      inc2 = 1
!
!   Initialisations for contributions from the Dirac-Coulomb
!   operator
!
      kt  = 0
      ipt = 1
!
      incor = 1

      ncoec = 0
!
      nctec   = 0

      if (ltrans) then

!        ...Initialisations for transverse interaction correction
         do 2 i = 1, orbital_count
            icore(i) = 0
            do j = 1, csf_count
               if (ichop (i,j) <= 0) goto 2
            end do
            icore(i) = 1
    2    continue

         nmcbp = 0
         ncore = 0
      end if

! Loop over rows of the Hamiltonian matrix - distributed

      do 10 ic = icstrt, csf_count, nprocs

         nelc = 0    ! counter - Number of non-zeros of this row

! Loop over columns of the current row

         irstart = 1
         do 85 ir = irstart, ic

! PER
            if (lfordr .and. (ir > iccut(1))) then
               if (ir/=ic) cycle
            end if
! PER

            elemnt = 0.0_real64     ! accumulates various contributions to H

!
!   Generate the integral list for the matrix element of the
!   one-body operators
!
            if (iprerun == 1) then
               inc1 = 0
               inc2 = 0
               if (ic<=ncsfpre .or. ic==ir) then
                  inc1 = 1
               end if
            end if

            if (iprerun == 2) then
!
!   Diagonal elements are always included
!   Off diagonal elements are included only if the
!   products of the weights from the prerun are larger
!   than the cutoff.
!
               if (ic == ir) then
                  inc1 = 1
                  inc2 = 1
               else
                  inc1 = 0
                  inc2 = 0
               end if
               do ipi = 1,nvec
                  precoeff =                                           &
                   dabs(evec1(ic+(ipi-1)*csf_count)*evec1(ir+(ipi-1)*csf_count))
                  if (precoeff > coeffcut1) inc1 = 1
                  if (precoeff > coeffcut2) inc2 = 1
               end do
            end if

!            ...INC1.EQ.1 ------------>
         if (inc1 == 1) then   !inc1 is always 1 without PRE-RUN
           call onescalar(ic,ir,ia,ib,tshell)
!
!   Accumulate the contribution from the one-body operators:
!   kinetic energy, electron-nucleus interaction; update the
!   angular integral counter
!
         if (ia /= 0) then
            if (ia == ib) then
               do ia = 1,orbital_count
                  iia = ia
                  tcoeff = tshell(ia)
                  if (abs (tcoeff) > cutoff) then
                     ncoec = ncoec + 1
                     call iabint (iia, iia, tegral)
                        !------------------------
                     elemnt = elemnt + tegral*tcoeff
                     if (lnms) then
                        call keint (iia,iia,tegral)
                        !------------------------
                        elemnt = elemnt + tegral*atwinv*tcoeff
                     end if
                     if (lvp) then
                        call vpint (iia, iia, tegral)
                        !------------------------
                        elemnt = elemnt + tegral*tcoeff
                     end if
                  end if
               end do
            else
               tcoeff = dble(tshell(1))
               if (dabs (tcoeff) > cutoff) then
                  ncoec = ncoec + 1
                  call iabint (ia, ib, tegral)
                        !------------------------
                  elemnt = elemnt + tegral*tcoeff
                  if (lnms) then
                     call keint (ia, ib, tegral)
                        !------------------------
                     elemnt = elemnt + tegral*atwinv*tcoeff
                  end if
                  if (lvp) then
                     call vpint (ia, ib, tegral)
                        !------------------------
                     elemnt = elemnt + tegral*tcoeff
                  end if
               end if
            end if
         end if
!
         ibug1 = 0
!
!   Accumulate the contributions from the two-electron
!   Coulomb operator and the mass polarisation; the latter
!   is computed first because the orbital indices may be
!   permuted by RKINTC
!
         angular_coefficient_count = 0
!
         call rkco_gg (ic, ir, cord, incor, 1)
!
         do 7 i = 1, angular_coefficient_count
            vcoeff = dble(angular_coefficients(i))
            if (dabs (vcoeff) > cutoff) then
               nctec = nctec + 1
               if (lsms) then
                  if (angular_coefficient_labels(5,i) == 1) then
                     call vint (angular_coefficient_labels(1,i), angular_coefficient_labels(3,i), tgrl1)
                     call vint (angular_coefficient_labels(2,i), angular_coefficient_labels(4,i), tgrl2)
                     elemnt = elemnt - tgrl1*tgrl2*atwinv*vcoeff
                  end if
               end if
               call rkintc (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),                  &
                            angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),                  &
                            angular_coefficient_labels(5,i), tegral)
               elemnt = elemnt + tegral*vcoeff
            end if
    7    continue
!
         ibug1 = 0

         end if  !inc1 is always 1 without PRE-RUN
!            ...INC1.EQ.1 <------------
!***********************************************************************
!            ...LTRANS .AND. (INC2.EQ.1) ------------>
         if (ltrans .and. (inc2==1)) then
            !IF (INC2 .EQ. 1) THEN  !inc2 is always 1 without PRE-RUN
!
!   Accumulate the contribution from the two-electron
!   transverse interaction operator
!
            angular_coefficient_count = 0
!
            call rkco_gg (ic, ir, breid, 1, 2)
!
            do 8 i = 1, angular_coefficient_count
               if (dabs (angular_coefficients(i)) > cutoff) then
                  nmcbp = nmcbp + 1
                  itype = abs (angular_coefficient_labels(6,i))
                  if (itype == 1) then
                     call brint1 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),               &
                                  angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                  angular_coefficient_labels(5,i), tegral)
                  else if (itype == 2) then
                     call brint2 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),               &
                                  angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                  angular_coefficient_labels(5,i), tegral)
                  else if (itype == 3) then
                     call brint3 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),               &
                                  angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                  angular_coefficient_labels(5,i), tegral)
                  else if (itype == 4) then
                     call brint4 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),               &
                                  angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                  angular_coefficient_labels(5,i), tegral)
                  else if (itype == 5) then
                     call brint5 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),               &
                                  angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                  angular_coefficient_labels(5,i), tegral)
                  else if (itype == 6) then
                     call brint6 (angular_coefficient_labels(1,i), angular_coefficient_labels(2,i),               &
                                  angular_coefficient_labels(3,i), angular_coefficient_labels(4,i),               &
                                  angular_coefficient_labels(5,i), tegral)
                  end if
                  contr = angular_coefficients(i)*tegral
                  if (angular_coefficient_labels(6,i) > 0) then
                     elemnt = elemnt + contr
                  else
!                        ...It comes here only when ic=ir=1
!                           clue: rkco<-breid<-talk<-label(6,i)
                     ncore = ncore + 1
                     elsto = elsto + contr
                  end if
               end if
    8       continue
!
            ibug1 = 0
!
!               ...ELSTO is a constant over all diagonals, thus its
!                  contribution to the total energy can be added later
!            IF (IR .EQ. IC) ELEMNT = ELEMNT + ELSTO
!
            !ENDIF   !inc2 is always 1 without PRE-RUN
         end if
!            ...LTRANS .AND. (INC2.EQ.1) <------------
!***********************************************************************
!
! Store this element if it is diagonal or its value is greater than
! CUTOFF
!
         if ( (ir == ic) .or. (dabs (elemnt) > cutoff) ) then
            nelc       = nelc + 1
            hamiltonian_values(nelc)  = elemnt
            matrix_partner_indices(nelc) = ir
         end if
!
   85    continue
! zou
!        print *, ic,SLF_EN(IC)
         if(lse) hamiltonian_values(nelc) = hamiltonian_values(nelc) + slf_en(ic)
! zou
!
!   This column is done; write it to disk
!
         write (imcdf) nelc, elsto, (hamiltonian_values(ir), ir = 1, nelc),       &
                                   (matrix_partner_indices(ir), ir = 1, nelc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This EAV (and the above EMT) does not have ELSTO.
         hamiltonian_energy_shift = hamiltonian_energy_shift + hamiltonian_values(nelc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   Output IC on screen to show how far the calculation has proceeded
!
!-----------------------------------------------------------------------
!cjb      if (mod(ic-1,100).eq.0) then
!cjb         PRINT '(A4,I10,A9,I10,A9,I5,A8,I5)', 'Row ', ic,       &
!cjb                 ' nnonz = ', nelc,                             &
!cjb                 ' block = ', jblock,' myid = ', myid
!cjb      ENDIF
!
!cjb          CALL CONVRT (IC,CNUM,LCNUM)
!cjb          if (mod(IC,100).eq.0) then
!cjb             PRINT *, 'Column '//CNUM(1:LCNUM)//' complete;'
!cjb          end if
!cjb
!
!cjb MPI progress begin --------------------------------------------------
!
!cjb just started
       if (ic <= nprocs) then
            print "(A6,I10,A9,I10,A9,I5,A8,I5)", "Start ", ic,       &
                    " nnonz = ", nelc,                             &
                    " block = ", jblock," myid = ", myid
        flush output_unit
        flush error_unit
       end if
!
!cjb progress
       if (ic > nprocs .and. ic < csf_count-nprocs .and. &
           mod (ic-1,100*nprocs) == myid) then
         if (myid == 0) then
            print "(A5,I10,A9,I10,A9,I5,A8,I5)", "Done ", ic,       &
                    " nnonz = ", nelc,                             &
                    " block = ", jblock," myid = ", myid
        flush output_unit
        flush error_unit
         else
            print "(A5,I10,A9,I10,A9,I5,A8,I5)", "Row  ", ic,       &
                    " nnonz = ", nelc,                             &
                    " block = ", jblock," myid = ", myid
        flush output_unit
        flush error_unit
         end if
       end if
!
!cjb almost finished
       if (ic > csf_count-nprocs) then

            print "(A7,I10,A9,I10,A9,I5,A8,I5)", "Finish ", ic,       &
                    " nnonz = ", nelc,                             &
                    " block = ", jblock," myid = ", myid
        flush output_unit
        flush error_unit
       end if
!cjb MPI progress end ----------------------------------------------------
!
!
!   Update the counter for the total number of elements
!
         stored_matrix_element_count = stored_matrix_element_count + nelc
!
   10 continue
!***********************************************************************
!
!   Deallocate storage for the arrays in /BUFFER/
!
      call alcbuf (3)

!     ...Locals
      call dalloc (hamiltonian_values, "EMT", "SETHAM")
      call dalloc (matrix_partner_indices, "IROW", "SETHAM")

!  Fill the common block /setham_to_genmat2/ for use in genmat2

      cutofftmp = cutoff
      ncoeitmp = ncoei
      ncoectmp = ncoec
      ncteitmp = rk_local_count(size(valteirk), myid, nprocs)
      nctectmp = nctec
      ntpitmp = ntpi
      nmcbptmp = nmcbp
      ncoretmp = ncore
      nvpitmp = nvpi
      nkeitmp = nkei
      nvintitmp = nvinti
      nelmnttmp = stored_matrix_element_count
      ncftmp = csf_count

      return
      end subroutine setham
