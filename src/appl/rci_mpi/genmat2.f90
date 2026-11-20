!***********************************************************************
!                                                                      *
      subroutine genmat2(irestart, nelmnt_a, elsto)
!
!   Get eav and do writings to the summary file .csum
!   The mpi version (genmat2mpi) also gets nelmnt_a and elsto
!
! Xinghong He 98-06-15
!
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use sparse_hamiltonian_storage, only: cutofftmp, ncftmp, ncoectmp, ncoeitmp, ncoretmp, nctectmp, ncteitmp, nelmnttmp, &
         & nkeitmp, nmcbptmp, ntpitmp, nvintitmp, nvpitmp  !! rather than setham_to_genmat2
      use hamiltonian_correction_flags, only: lfordr, lnms, lse, lsms, ltrans, lvp
      use eigenstate_expansions, only: hamiltonian_energy_shift, eavff, eavii, shifted_eigenvalues, evalff, evalii, &
          csf_mixing_coefficients, evecff, evecii
      use mpi_state, only: mpi_allreduce, mpi_bcast, mpi_comm_world, mpi_double_precision, mpi_error_code, &
         & mpi_integer, mpi_integer8, mpi_max, mpi_rank, mpi_reduce, mpi_sum
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)        :: irestart
      integer(int64), intent(out) :: nelmnt_a
      real(real64), intent(in)   :: elsto
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ntpi_a(6)
      integer               :: ncoei_a, nctei_a, nctec_a, nmcbp_a,    &
                               ncore_a, nvpi_a, nkei_a, nvinti_a,     &
                               ncoec_a
      real(real64)          :: densty, eav_a
!-----------------------------------------------
!-----------------------------------------------------------------------
! Transfer parameters around and print out on node-0.
!  .Temporary variables (ending with "_a") are used to store the
!   sums over nodes;
!  .Outputs to stream 24 are done on node-0;
! The following parameters are accumulated in setham which will not
! contain the correct values in restart mode. And thus is skipped.
!-----------------------------------------------------------------------
      if (irestart /= 0) then                    ! non-restart mode
!         WRITE (24, 301) CUTOFFTMP
!         WRITE (24, 302) NCOEITMP
!         WRITE (24, 303) NCOECTMP
!         WRITE (24, 304) NCTEITMP
!         WRITE (24, 305) NCTECTMP
!         IF (LTRANS) THEN
!            WRITE (24, 306) NTPITMP
!            WRITE (24, 307) NMCBPTMP
!            WRITE (24, 308) NCORETMP
!         ENDIF
!         IF (LVP) WRITE (24, 309) NVPITMP
!         IF (LNMS) WRITE (24, 310) NKEITMP
!         IF (LSMS) WRITE (24, 311) NVINTITMP
!      ELSE
!         WRITE (24, *) 'Restart mode --- no report on radial integrals'
      call mpi_reduce (ncoeitmp, ncoei_a, 1, mpi_integer, mpi_sum, 0,  &
                                       mpi_comm_world, mpi_error_code)
      call mpi_reduce (ncoectmp, ncoec_a, 1, mpi_integer, mpi_sum, 0,  &
                                       mpi_comm_world, mpi_error_code)
      call mpi_reduce (ncteitmp, nctei_a, 1, mpi_integer, mpi_sum, 0,  &
                                       mpi_comm_world, mpi_error_code)
      call mpi_reduce (nctectmp, nctec_a, 1, mpi_integer, mpi_sum, 0,  &
                                       mpi_comm_world, mpi_error_code)

      if (ltrans) then
         call mpi_reduce (ntpitmp, ntpi_a, 6, mpi_integer, mpi_sum, 0, &
                                       mpi_comm_world, mpi_error_code)
         call mpi_reduce (nmcbptmp, nmcbp_a, 1, mpi_integer, mpi_sum,0,&
                                       mpi_comm_world, mpi_error_code)
         call mpi_reduce (ncoretmp, ncore_a, 1, mpi_integer,mpi_sum,0, &
                                       mpi_comm_world, mpi_error_code)
      end if

      if (lvp)                                                         &
         call mpi_reduce (nvpitmp, nvpi_a, 1, mpi_integer, mpi_sum, 0, &
                                       mpi_comm_world, mpi_error_code)

      if (lnms)                                                        &
         call mpi_reduce (nkeitmp, nkei_a, 1, mpi_integer, mpi_sum, 0, &
                                       mpi_comm_world, mpi_error_code)

      if (lsms)                                                        &
         call mpi_reduce (nvintitmp, nvinti_a, 1, mpi_integer, mpi_sum,&
                                       0, mpi_comm_world, mpi_error_code)

      if (mpi_rank == 0) then
         write (24,'("CUTOFF set to ",1p,d17.10)') cutofftmp
         write (24,'("Dirac-Coulomb one-e radial integrals:",1i8)') ncoei_a
         write (24,'("One-e angular integrals that exceed CUTOFF: ",1i8)') ncoec_a
         write (24,'("Coulomb two-e radial integrals: ",1i8)') nctei_a
         write (24,'("Two-e angular integrals that exceed CUTOFF: ",1i11)') nctec_a
         if (ltrans) then
            write (24,'("Transverse two-e radial integrals: "/,6i8)') ntpi_a
            write (24,'("MCBP coefficients that exceed CUTOFF: ",1i16)') nmcbp_a
            write (24,'("Core coefficients that exceed CUTOFF: ",1i8)') ncore_a
         end if
         if (lvp) write (24,'("Vacuum polarisation integrals: ",1i8)') nvpi_a
         if (lnms) write (24,'("Kinetic energy integrals: ",1i8)') nkei_a
         if (lsms) write (24,'("Vinti integrals: ",1i8)') nvinti_a
      end if
            else
      if (mpi_rank == 0) then
         write (24,*) "Restart mode --- no report on radial integrals"
      end if
      end if                                      !(irestart .NE. 0) ! non-restart mode



!-----------------------------------------------------------------------
! ELSTO, EAV are not only for print-out, but also used later.
! density of the Hamiltonian matrix is only for print-out.
! At this place Hamiltonian matrix elements (EMT on each node) do
!   _not_ contain ELSTO. ELSTO will be added to the total energy
!  later with EAV.
!-----------------------------------------------------------------------

      call mpi_allreduce (nelmnttmp, nelmnt_a, 1, mpi_integer8,  &
!AK:      MPI_INTEGER8 required to match the sizes of NELMNT_a and NELMNTtmp
!CFF &    MPI_SUM, MPI_COMM_WORLD, ierr) ! want MAX not SUM
     &    mpi_max, mpi_comm_world, mpi_error_code)
!cjb  NELMNT_a decides between sparse and dense
!cjb  see comments in maneig.f90

      call mpi_bcast (elsto, 1, mpi_double_precision, 0,         &
                            mpi_comm_world, mpi_error_code)
      call mpi_allreduce (hamiltonian_energy_shift, eav_a, 1, mpi_double_precision,   &
                              mpi_sum, mpi_comm_world, mpi_error_code)
      hamiltonian_energy_shift = eav_a / dble (ncftmp) + elsto

      if (mpi_rank == 0) then
         write (24,'("Elements that exceed CUTOFF in the lower"," triangle of the H matrix: ",1i11)') nelmnt_a
         write (24, *)
         write (24, '("Average energy = ",1p,d19.12," Hartrees.")') hamiltonian_energy_shift
      end if

!cjb 1I8 -> 1I16
! 307 FORMAT('MCBP coefficients that exceed CUTOFF: ',1I8)

      return
      end subroutine genmat2
