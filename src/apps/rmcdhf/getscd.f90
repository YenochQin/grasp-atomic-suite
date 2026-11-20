!***********************************************************************
!                                                                      *
      subroutine getscd(extended_optimal_level, block_labels, isotope_file_path, radial_file_path, nuclear_polarization_file_path)
      use, intrinsic :: iso_fortran_env, only: input_unit
      use calculation_file_units, only: calculation_log_unit
!                                                                      *
!   Interactively determines the data governing the SCF problem.       *
!                                                                      *
!   Call(s) to: [LIB92]: CONVRT, GETRSL, GETYN, NUCPOT, RADGRD,        *
!                        SETISO, SETQIC, SETRWF.                       *
!               [RSCF92]: GETALD, GETOLD
!                                                                      *
!   Written by Farid A. Parpia            Last revision: 27 Dec 1992   *
!
! Xinghong He 98-08-06
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,    only: max_orbitals
      use radial_node_count_settings, only: count_context, count_trace_rank, thresh
      use scf_damping_factors, only: orbital_damping_factors, mixing_damping_factors
      use atomic_calculation_state, only: accy, max_radial_points, ncmin, nscf, nsic, nsolv, nuclear_charge, speed_of_light_au, &
          vacuum_speed_of_light_au
      use calculation_input_defaults, only: inpci, ndef, ndump
      use fixed_orbital_flags, only: is_orbital_fixed
      use standard_io_units, only: stderr_unit, stdin_unit, stdout_unit
      use radial_grid_state, only: transformed_grid_step, asymptotic_grid_step, radial_point_count, radial_grid_scale, radial_grid
      use orbital_inversion_flags, only: noinvt
      use rmcdhf_mcp_options, only: diag, lfordr
      use mpi_s, only: mpi_rank
      use expected_orbital_nodes, only: nnodep
      use nuclear_model_parameters_state, only: nuclear_model_parameter_count, nuclear_model_parameters
      use nuclear_potential_state, only: negative_radius_times_potential
      use nuclear_polarization_state, only: npol_message_len
      use nuclear_polarization_initialization, only: npol_init
      use orbital_csf_basis , only: orbital_kappa, orbital_symmetry_labels, orbital_principal_n, orbital_count, ped
      use scf_orthogonality_control, only: orthst
      use scf_iteration_state, only: scnsty, method
      use radial_wavefunctions, only: max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
          pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use yes_no_input, only: getyn
      use isotope_data_input, only: setiso
      use radial_stencil_weights, only: setqic
      use radial_grid_setup, only: radgrd
      use nuclear_coulomb_potential, only: nucpot
      use radial_wavefunction_input, only: setrwfa
      use getald_i, only: getald
      use getold_i, only: getold
      use integer_format, only: convrt
      use subshell_label_selection, only: getrsl
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical, intent(out) :: extended_optimal_level
      character  :: isotope_file_path*(*)
      character  :: radial_file_path*(*)
      character(*), intent(in) :: nuclear_polarization_file_path
      character  :: block_labels(*)*8
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      character(len=6), parameter :: routine_name = "GETSCD"
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: np_status
      character(npol_message_len) :: np_message
      integer :: iqadum
      integer  :: indx(max_orbitals)
      integer :: i, j, iend, ibeg, text_length, nsubs, loc
      real(real64) :: odampu, cdampu
      logical :: yes
      character :: record*80, number_text*20
!-----------------------------------------------
!
!   Open, check, load data from, and close the  .iso  file
!
      call setiso (isotope_file_path)
!
!   Set default speed of light and grid parameters
!
      speed_of_light_au = vacuum_speed_of_light_au
      if (nuclear_model_parameter_count == 0) then
         radial_grid_scale = exp((-65.0_real64/16.0_real64))/nuclear_charge
         transformed_grid_step = 0.5_real64**4
         radial_point_count = min(220,max_radial_points)
      else
!         default comes here
!CFF     .. should be Z-dependent
         radial_grid_scale = 2.0e-6_real64/nuclear_charge
         transformed_grid_step = 5.0e-2_real64
         radial_point_count = max_radial_points
      end if
      asymptotic_grid_step = 0.0_real64
!
!   ACCY is an estimate of the accuracy of the numerical procedures
!
      accy = transformed_grid_step**6
      if (ndef /= 0) then

         write (stderr_unit,"(A)",advance="NO")"Change the default speed of",&
            " light or radial grid parameters?  (y/n) "
         yes = getyn()
         if (yes) then
            write (stderr_unit,*) "Speed of light = ",vacuum_speed_of_light_au,";", " revise ?"
            yes = getyn()
            if (yes) then
               write (stderr_unit, *) "Enter the revised value:"
               read (input_unit, *) speed_of_light_au
            end if
!
!   Determine the parameters controlling the radial grid
!
            write (stderr_unit, *) "The default radial grid parameters for ",&
               "this case are:"
            write (stderr_unit, *) " RNT = ", radial_grid_scale
            write (stderr_unit, *) " H = ", transformed_grid_step
            write (stderr_unit, *) " HP = ", asymptotic_grid_step
            write (stderr_unit, *) " N = ", radial_point_count
            write (stderr_unit, *) " revise these values?"
            yes = getyn()
            if (yes) then
               write (stderr_unit, *) "Enter RNT:"
               read (input_unit, *) radial_grid_scale
               write (stderr_unit, *) "Enter H:"
               read (input_unit, *) transformed_grid_step
               write (stderr_unit, *) "Enter HP:"
               read (input_unit, *) asymptotic_grid_step
               write (stderr_unit, *) "Enter N:"
               read (input_unit, *) radial_point_count
!b Revised grid
               write (stderr_unit,*) "Revised RNT = ", radial_grid_scale
               write (stderr_unit,*) "Revised H   = ", transformed_grid_step
               write (stderr_unit,*) "Revised HP  = ", asymptotic_grid_step
               write (stderr_unit,*) "Revised N   = ", radial_point_count
            end if
         end if
!
!b
!b read ACCY on input
!b
         write (stderr_unit,*) "Revise the default ACCY = ", accy
         yes = getyn ()
            if (yes) then
               write (stderr_unit,*) "Enter ACCY:"
               read *, accy
!b Revised ACCY
               write (stderr_unit,*) "Revised ACCY = ", accy
         end if
      end if
!
!   Set up the coefficients for the numerical procedures
!
      call setqic
!
!   Generate the radial grid and all associated arrays
!
      call radgrd
!
!   Generate $- r \times V_ (r)$
!
      call nucpot
      call npol_init(nuclear_polarization_file_path, nuclear_charge, radial_point_count, radial_grid, &
         & negative_radius_times_potential, &
         & np_status, np_message)
      if (np_status /= 0) then
         write(stderr_unit,"(A)") trim(np_message)
         error stop 1
      end if
!
!   Load the subshell radial wavefunction estimates
!
      call setrwfa (radial_file_path)
!
!   Set some defaults
!
      thresh = 0.05_real64

!      IORDER(I) = I    ! Completely determined in GETOLD
      method(:orbital_count) = 3
      noinvt(:orbital_count) = .true.
!CFF      ODAMP(:NW) = 1.D0
      orbital_damping_factors(:orbital_count) = 0.0_real64
      ped(:orbital_count)   = 0.0_real64
      scnsty(:orbital_count) = 0.0_real64

      where (orbital_kappa(:orbital_count) < 0)
         nnodep(:orbital_count) = orbital_principal_n(:orbital_count) + orbital_kappa(:orbital_count)
      else where
         nnodep(:orbital_count) = orbital_principal_n(:orbital_count) - orbital_kappa(:orbital_count) - 1
      end where

      if (diag) then
         extended_optimal_level = .false.
         call getald                             ! (E)AL type calculation,
                              ! H(DC) will not be diagonalised
      else if (lfordr) then
         extended_optimal_level = .true.
         call getold (block_labels)                     ! EOL type calculation
      else
!GG         WRITE (ISTDE, '(A)', ADVANCE='NO') '(E)OL type calculation?  (y/n) '
         extended_optimal_level = .true.
!GG         EOL = GETYN()
         if (extended_optimal_level) then
            call getold (block_labels)
         else
            call getald                          ! (E)AL type calculation,
                              ! H(DC) will not be diagonalised
         end if
      end if

      write (stderr_unit, *) "Enter the maximum number of SCF cycles:"
      read (*, *) nscf

      if (ndef==0) then
         write(calculation_log_unit,*) nscf,"! Number of SCF cycles"
      end if
!
!   Allow the user to modify other defaults
!
      if (ndef /= 0) then
         write (stderr_unit, "(A)", advance="NO") "Modify other defaults?  (y/n) "
         yes = getyn()
      else
         yes = .false.
      end if

      if (.not.yes) return
!=======================================================================
! From here to end, "other defaults" are handled. For simplicity
! We'll let node-0 do the job and then broadcast results to all
! nodes.
!=======================================================================
                  !-------------------------------------------
      if (mpi_rank == 0) then                        ! This is a _big_ IF
                  !-------------------------------------------
!
!   THRESH
!
         write (stderr_unit,*)"An oscillation in the large-component of the ",&
            "radial wavefunction is diregarded"
         write (stderr_unit, *) "for the purposes of node counting if its ", &
            "amplitude is less than 1/20 the"
         write (stderr_unit, *) "maximum amplitude.   Revise this?"
         yes = getyn()
         if (yes) then
    3       continue
            write (stderr_unit, *) "Enter the new threshold value:"
            read (*, *) thresh
            if (thresh <= 0.0_real64) then
               write (stderr_unit, *) routine_name, ": This must exceed 0;"
               goto 3
            end if
         end if
         yes = .false.
!
!   METHOD
!

! Piece only for printing...

         do i = 1, 4

            do j = 1, orbital_count
               if (.not.(method(j)==i .and. .not.is_orbital_fixed(j))) cycle
               write (stderr_unit, *) "Method ", i, " is used for ", &
                  "integrating the radial differential ", &
                  "equation for subshells"
               goto 9
            end do

            cycle

    9       continue
            iend = 0
            do j = 1, orbital_count
               if (method(j)==i .and. .not.is_orbital_fixed(j)) then
                  ibeg = iend + 1
                  iend = ibeg
                  record(ibeg:iend) = " "
                  call convrt (orbital_principal_n(j), number_text, text_length)
                  ibeg = iend + 1
                  iend = ibeg + text_length - 1
                  record(ibeg:iend) = number_text(1:text_length)
                  ibeg = iend + 1
                  if (orbital_kappa(j) < 0) then
                     iend = ibeg
                     record(ibeg:iend) = orbital_symmetry_labels(j)(1:1)
                  else
                     iend = ibeg + 1
                     record(ibeg:iend) = orbital_symmetry_labels(j)(1:2)
                  end if
               end if
               if (iend <= 76) cycle
               write (stderr_unit, *) record(1:iend)
               iend = 0
            end do
            if (iend<=0 .or. mpi_rank/=0) cycle
            write (stderr_unit, *) record(1:iend)
         end do

! Reads user inputs and fills array indx(1:nsubs) where nsubs itself
! is an output from getrsl.
! METHOD(1:4) is the only output. indx() and nsubs are discarded

         write (stderr_unit, *) "Select a different integration method for ", &
            "any subshell radial wavefunction?"
         yes = getyn()
         if (yes) then
            do i = 1, 4
               write (stderr_unit, *) "Method ", i, ":"
               call getrsl (indx, nsubs)
               do j = 1, nsubs
                  loc = indx(j)
                  if (is_orbital_fixed(loc)) cycle
                  method(loc) = i
               end do
            end do
         end if
!
!   NOINVT
!
         write (stderr_unit, *) "The first oscillation of the large component"

         do i = 1, orbital_count
            if (.not.(noinvt(i) .and. .not.is_orbital_fixed(i))) cycle
            write (stderr_unit, *) "of the following radial wavefunctions ", &
               "will be required to be positive"
            goto 15
         end do

         write (stderr_unit, *) "of all radial wavefunctions will be required ", &
            "to be positive.   Revise this?"
         yes = getyn()
         goto 17
   15    continue
         iend = 0
         do i = 1, orbital_count
            if (noinvt(i) .and. .not.is_orbital_fixed(i)) then
               ibeg = iend + 1
               iend = ibeg
               record(ibeg:iend) = " "
               call convrt (orbital_principal_n(i), number_text, text_length)
               ibeg = iend + 1
               iend = ibeg + text_length - 1
               record(ibeg:iend) = number_text(1:text_length)
               ibeg = iend + 1
               if (orbital_kappa(i) < 0) then
                  iend = ibeg
                  record(ibeg:iend) = orbital_symmetry_labels(i)(1:1)
               else
                  iend = ibeg + 1
                  record(ibeg:iend) = orbital_symmetry_labels(i)(1:2)
               end if
            end if
            if (iend <= 76) cycle
            write (stderr_unit, *) record(1:iend)
            iend = 0
         end do
         if (iend > 0) write (stderr_unit, *) record(1:iend)
         write (stderr_unit, *) "Revise this?"
         yes = getyn()
   17    continue
         if (yes) then
            write (stderr_unit, *) "Suppressing enforcement of positive first ", &
               "oscillation:"
            call getrsl (indx, nsubs)
            do i = 1, nsubs
               loc = indx(i)
               if (is_orbital_fixed(loc)) cycle
               noinvt(loc) = .true.
            end do
         end if
!
!   ODAMP
!
         do i = 1, orbital_count
            if (.not.(orbital_damping_factors(i)/=0.0_real64 .and. .not.is_orbital_fixed(i))) cycle
            write (stderr_unit, *) "Subshell accelerating parameters have ", &
               "been set.   Revise these?"
            yes = getyn()
            goto 20
         end do
         write (stderr_unit, *) "Set accelerating parameters for subshell ", &
            "radial wavefunctions?"
         yes = getyn()
   20    continue
         if (yes) then
            write (stderr_unit, *) "Different accelerating parameters for ", &
               "different subshell radial wavefunction?"
            yes = getyn()
            if (yes) then
   21          continue
               write (stderr_unit, *) "Enter an accelerating parameter"
               write (stderr_unit, *) " (0< ODAMP < 1 allows ODAMP to be ", &
                  "reduced as convergence is approached;"
               write (stderr_unit, *) " -1 < ODAMP < 0 implies |ODAMP| is ", &
                  "held constant):"
               read (*, *) odampu
               if (abs(odampu)==0.0_real64 .or. abs(odampu)>=1.0_real64) then
                  write (stderr_unit, *) routine_name, ": Value out of range ..."
                  goto 21
               else
                  call getrsl (indx, nsubs)
                  do i = 1, nsubs
                     loc = indx(i)
                     if (is_orbital_fixed(loc)) cycle
                     orbital_damping_factors(loc) = odampu
                  end do
               end if
            else
   23          continue
               write (stderr_unit, *) "Enter the accelerating parameter"
               write (stderr_unit, *) " (0< ODAMP < 1 allows ODAMP to be ", &
                  "reduced as convergence is approached;"
               write (stderr_unit, *) " -1 < ODAMP < 0 implies |ODAMP| is ", &
                  "held constant):"
               read (*, *) odampu
               if (abs(odampu)==0.0_real64 .or. abs(odampu)>=1.0_real64) then
                  write (stderr_unit, *) routine_name, ": Value out of range ..."
                  goto 23
               else
                  where (.not.is_orbital_fixed(:orbital_count))
                     orbital_damping_factors(:orbital_count) = odampu
                  end where
               end if
            end if
         end if
!
!   CDAMP
!
         write (stderr_unit, *) "Set accelerating parameters for the ", &
            "eigenvectors?"
         yes = getyn()
         if (yes) then
            write (stderr_unit, *) "Different accelerating parameters ", &
               "for each eigenvector?"
            yes = getyn()
            if (yes) then
               write (stderr_unit, *) "Enter an accelerating parameter for"
               call convrt (ncmin, record, text_length)
               write (stderr_unit, *) " each of the "//record(1:text_length)//" levels :"
               read (*, *) (mixing_damping_factors(i),i=1,ncmin)
            else
               write (stderr_unit, *) "Enter the accelerating parameter:"
               read (*, *) cdampu
               mixing_damping_factors(:ncmin) = cdampu
            end if
         end if
!
!   NSIC
!
         write (stderr_unit, *) "Following the improvement of each of the ", &
            "subshell radial wavefunctions in turn, "
         write (stderr_unit, *) "the ", nsic, " least self-consistent", &
            " functions will be improved at the"
         write (stderr_unit, *) "end of the first SCF cycle. Revise this ", &
            "setting?"
         yes = getyn()
         if (yes) then
            write (stderr_unit, *) "Enter the number of additional ", "improvements:"
            read (*, *) nsic
         end if
!
!   NSOLV
!
         write (stderr_unit, *) "The maximum number of cycles in attempting ", &
            "to solve each radial equation is "
         write (stderr_unit, *) nsolv, " times the principal quantum", &
            " number of the radial"
         write (stderr_unit, *) "wave-function to be estimated.   ", &
            "Revise this setting?"
         yes = getyn()
         if (yes) then
            write (stderr_unit, *) "Enter the factor that multiplies the ", &
               "principal quantum number:"
            read (*, *) nsolv
         end if
!
!   Orthogonalisation
!
         if (orthst) then
            write (stderr_unit, *) "Subshell radial wavefunctions will be ", &
               "Schmidt orthogonalised immediately"
            write (stderr_unit, *) "following their estimation to all ", &
               "functions with poorer self-consistency."
            write (stderr_unit, *) " Revise this?"
            yes = getyn()
            if (yes) orthst = .false.
         else
            write (stderr_unit, *) "Subshell radial wavefunctions will be ", &
               "Schmidt orthogonalised at the end of"
            write (stderr_unit, *) "each SCF cycle.   Revise this?"
            yes = getyn()
            if (yes) orthst = .true.
         end if
                  !-------------------------------------------
      end if                                      ! end of the _big_ IF
                  !-------------------------------------------
      return
      end subroutine getscd
