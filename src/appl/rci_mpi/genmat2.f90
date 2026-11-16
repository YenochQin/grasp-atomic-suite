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
!...Translated by Pacific-Sierra Research 77to90  4.3E  14:04:58   1/ 3/07
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64, int64
      use hmat_c  !! rather than setham_to_genmat2
      use decide_c
      use eigv_c
      use mpi_c
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
      endif

      if (lvp)                                                         &
         call mpi_reduce (nvpitmp, nvpi_a, 1, mpi_integer, mpi_sum, 0, &
                                       mpi_comm_world, mpi_error_code)

      if (lnms)                                                        &
         call mpi_reduce (nkeitmp, nkei_a, 1, mpi_integer, mpi_sum, 0, &
                                       mpi_comm_world, mpi_error_code)

      if (lsms)                                                        &
         call mpi_reduce (nvintitmp, nvinti_a, 1, mpi_integer, mpi_sum,&
                                       0, mpi_comm_world, mpi_error_code)

      if (mpi_rank .eq. 0) then
         write (24,301) cutofftmp
         write (24,302) ncoei_a
         write (24,303) ncoec_a
         write (24,304) nctei_a
         write (24,305) nctec_a
         if (ltrans) then
            write (24,306) ntpi_a
            write (24,307) nmcbp_a
            write (24,308) ncore_a
         endif
         if (lvp) write (24,309) nvpi_a
         if (lnms) write (24,310) nkei_a
         if (lsms) write (24,311) nvinti_a
      endif
            else
      if (mpi_rank .eq. 0) then
         write (24,*) 'Restart mode --- no report on radial integrals'
      endif
      endif                                      !(irestart .NE. 0) ! non-restart mode



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

      if (mpi_rank .eq. 0) then
         write (24,312) nelmnt_a
         write (24, *)
         write (24, 300) hamiltonian_energy_shift
      endif

  300 format('Average energy = ',1p,d19.12,' Hartrees.')
  301 format('CUTOFF set to ',1p,d17.10)
  302 format('Dirac-Coulomb one-e radial integrals:',1i8)
  303 format('One-e angular integrals that exceed CUTOFF: ',1i8)
  304 format('Coulomb two-e radial integrals: ',1i8)
  305 format('Two-e angular integrals that exceed CUTOFF: ',1i11)
  306 format('Transverse two-e radial integrals: '/,6i8)
!cjb 1I8 -> 1I16
! 307 FORMAT('MCBP coefficients that exceed CUTOFF: ',1I8)
  307 format('MCBP coefficients that exceed CUTOFF: ',1i16)
  308 format('Core coefficients that exceed CUTOFF: ',1i8)
  309 format('Vacuum polarisation integrals: ',1i8)
  310 format('Kinetic energy integrals: ',1i8)
  311 format('Vinti integrals: ',1i8)
  312 format('Elements that exceed CUTOFF in the lower',&
         ' triangle of the H matrix: ',1i11)
  313 format('Density of the H(amiltonian) matrix: ',1p,d22.15)

      return
      end subroutine genmat2
