!***********************************************************************
      subroutine improvmpi (eol, j, lsort, dampmx)
!   The difference from the serial version is that it calls MPI        *
!   version subroutines (setlagmpi, cofpotmpi, matrixmpi, newcompi).   *
!                                                                      *
!   Improve the orbital J.                                             *
!                                                                      *
!   Call(s) to: [RSCF92]: DACON, DAMPCK, DAMPOR, LAGCON, matrixmpi,    *
!                         newcompi, ORTHOR, ROTATE, SETCOF, setlagmpi, *
!                         SOLVE, XPOT, YPOT.                           *
!               [LIB92]: ORTHSC, QUAD.                                 *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 22 Dec 1992   *
!   Modified by Xinghong He                 Last update: 05 Aug 1988   *
!   Modified for ifort -i8 by A. Kramida (AK) Last update 22 Mar 2016  *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use memory_man
      use corre_c
      use damp_c
      use def_c
      use grid_c
      use int_c
      use mpi_c
      use orb_c
      use orthct_c
      use scf_c
      use tatb_c
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use cofpotmpi_i
      use defcor_i
      use solve_i
      use orthsc_i
!GG      USE matrixmpi_I
!GG      USE newcompi_I
      use setlagmpi_i
      use quad_i
      use consis_i
      use dampck_i
      use dampor_i
      use orthy_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: j
      real(real64), intent(inout) :: dampmx
      logical, intent(in) :: eol
      logical :: lsort
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: p2 = 2.0e-1_real64
      real(real64), parameter :: p005 = 5.0e-3_real64
      real(real64), parameter :: p0001 = 1.0e-4_real64
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: ipr, npts, inv, jp, nnp, i, nwww
      integer :: i_mpi, ndcof_max, ntot, i_last, iproc, ndcip, jproc
      integer :: ifound, ind_buf, k
      real(real64) :: ed1, gamaj, ed2, wtaev, dnorm, dnfac, del1, del2, odampj
      logical :: fail, first
      real(real64), pointer :: da_buffer(:)
      integer, pointer :: nda_buffer(:),ndcof_buffer(:)
!-----------------------------------------------
!
!   C Froese Fischer's IPR and ED1 parameter
!
      data ipr/ 0/
      data ed1/ 0._real64/
      data first/ .false./
!
!
!-----------------------------------------------------------------------
!AK Handling the -i8 option of ifort and -fdefault-integer-8 option of gfortran
!      ISIZE = sizeof(NCF)
      i_mpi = mpi_integer
!      if (ISIZE.EQ.8) I_MPI = MPI_INTEGER8
!
      gamaj = orbital_origin_exponents(j)
!
!   C Froese Fischer's parameters IPR, ED1, ED2 are set and
!   used in this routine and in DAMPCK
!
    1 continue
      ed2 = orbital_energy_parameters(j)
!
!   Set up the exchange potential and arrays XU, XV as appropriate
!
!   Set coefficients for YPOT, XPOT, DACON
!   Compute direct potential, exchange potential
!   Add in Lagrange-multiplier contribution
!   Add in derivative-terms contribution
!
      npts = radial_point_count
      call cofpotmpi (eol, j, npts)
!
!   Calculate deferred corrections
!
      call defcor (j)
!
!   Solve the Dirac equation
!
!AK     NDCOF, NDA, and DA cannot be reduced by finding the max of NDCOF
!          and summing up DA from all nodes. The correct replacement is below
!cjb MPI_ALLREDUCE NDCOF
!      call MPI_ALLREDUCE(ndcof,ndcof_buffer,1,
!     :          I_MPI,MPI_MAX,MPI_COMM_WORLD,ierror)
!      if (ndcof_buffer .gt. ndcof) then
!        print *, 'improvmpi: ndcof, ndcof_buffer, myid=',
!     &                     ndcof, ndcof_buffer, myid
!        write (*,'(A35,4I10)') 'improv11: myid,NDCOF,ndcbfr=',
!     &    myid,NDCOF,ndcof_buffer
!        do indcof = ndcof+1, ndcof_buffer
!           da(indcof) = 0.0
!           nda(indcof) = 0
!        enddo
!        ndcof = ndcof_buffer
!      endif
      inv = 0
!cjb MPI_ALLREDUCE DA
!      if (ndcof.gt.0) then
!         call alloc(pda_buffer, ndcof, 8)
!         call alloc(pnda_buffer, ndcof, ISIZE)
!         call MPI_ALLREDUCE(da,da_buffer,ndcof,
!     :          MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
!         da(1:ndcof) = da_buffer(1:ndcof)
!      nda(1:ndcof) = nda_buffer(1:ndcof)
!         call dalloc(pda_buffer)
!         call dalloc(pnda_buffer)
!      ENDIF
!AK      Consolidate NDA and DA from all nodes
      if (mpi_process_count.gt.1) then
         call alloc(ndcof_buffer,mpi_process_count,'ndcof_buffer','IMPROVmpi')
         call mpi_gather(ndcof, 1, i_mpi, ndcof_buffer,            &
              1, i_mpi, 0, mpi_comm_world, mpi_error_code)
         call mpi_bcast (ndcof_buffer, mpi_process_count, i_mpi, 0,           &
                               mpi_comm_world, mpi_error_code)
         ndcof_max = 0
         do i = 1, mpi_process_count
            if (ndcof_buffer(i).gt.ndcof_max) ndcof_max=ndcof_buffer(i)
         enddo
         if (ndcof_max.gt.0) then
            ntot = ndcof_max*mpi_process_count
            call alloc(nda_buffer,ntot,'nda_buffer','IMPROVmpi')
            call mpi_gather(nda, ndcof_max, i_mpi, nda_buffer,     &
                 ndcof_max, i_mpi, 0, mpi_comm_world, mpi_error_code)
            call mpi_bcast (nda_buffer, ntot, i_mpi, 0,            &
                               mpi_comm_world, mpi_error_code)
            call alloc(da_buffer, ntot, 'da_buffer','IMPROVmpi')
            call mpi_gather(da, ndcof_max, mpi_double_precision,   &
              da_buffer, ndcof_max, mpi_double_precision, 0,       &
              mpi_comm_world, mpi_error_code)
            call mpi_bcast (da_buffer, ntot, mpi_double_precision, &
                               0, mpi_comm_world, mpi_error_code)
            i_last = 0
            do iproc = 1, mpi_process_count
               ndcip = ndcof_buffer(iproc)
               do jproc = 1, ndcip
                  ifound = 0
                  ind_buf = (iproc - 1)*ndcof_max + jproc
                  do k = 1, i_last
                     if (nda_buffer(ind_buf) .eq. nda(k)) then
                        da(k) = da(k) + da_buffer(ind_buf)
                        ifound = k
                        exit
                     endif
                  enddo
                  if (ifound.eq.0) then
                     i_last = i_last + 1
                     do while (i_last .gt. nddim)
                        if (nddim .gt. 0) then
                           nddim = 2*nddim
!GG                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 2)
                           call ralloc(nda,nddim,'NDA','IMPROVmpi')
                        else
                          nddim = 64
!GG                           CALL alcsca (PNTNDA, PNTRDA, NDDIM, 1)
                           call alloc(nda,nddim,'NDA','IMPROVmpi')
                        endif
                     enddo
                     nda(i_last) = nda_buffer(ind_buf)
                     da(i_last) = da_buffer(ind_buf)
                  endif
               enddo
            enddo
            ndcof = ndcof_max
            call dalloc(nda_buffer,'nda_buffer','IMPROVmpi')
            call dalloc(da_buffer,'da_buffer','IMPROVMPI')
         endif
         call dalloc(ndcof_buffer,'ndcof_buffer','IMPROVmpi')
      endif

      call solve (j, fail, inv, jp, nnp)
!
!   Upon failure issue message; take corrective action if possible
!
      if (fail) then
         if (mpi_rank == 0) write (*, 300) orbital_principal_n(j), orbital_symmetry_labels(j), method(j)
         if (method(j) /= 2) then
            method(j) = 2
!XHH orthsc does not have any argument
!    Orbital J [PF() and QF()]is not updated, why redo orthogonalization
            call orthsc
!CFF        ... avoid rediagonalization
!            IF (EOL) THEN
!               CALL MATRIXMPI
!               CALL NEWCOMPI (WTAEV)
!            ENDIF
            call setlagmpi (eol)
            go to 1
         else
            if (mpi_rank == 0) write (*, 301)
            !CALL TIMER (0)
            error stop
         endif
      endif
!
!   Compute norm of radial function
!
      ta(1) = 0._real64
      ta(2:mtp0) = (p(2:mtp0)**2+q(2:mtp0)**2)*radial_grid_derivative(2:mtp0)
      mtp = mtp0

      call quad (dnorm)

!   Determine self-consistency [multiplied by SQRT(UCF(J))]

      call consis (j)
!
!   Normalize
!
      dnfac = 1._real64/dsqrt(dnorm)
      p0 = p0*dnfac
      p(:mtp0) = p(:mtp0)*dnfac
      q(:mtp0) = q(:mtp0)*dnfac
!
!   Check if different method should be used or if improvement
!   count should be reduced
!
      del1 = dabs(1._real64 - ed2/orbital_energy_parameters(j))
      if (method(j) == 1) then
         del2 = dmax1(dabs(1._real64 - dsqrt(dnorm)),dabs(dnfac - 1._real64))
         if (del1<p005 .and. del2>p2) then
            method(j) = 2
            go to 1
         endif
      else
         if (del1<p0001 .and. nsic>1) nsic = nsic - 1
      endif
!
!   Damp the orbital --- if not converged
!
      if (scnsty(j) > accy) then
         call dampck (ipr, j, ed1, ed2)
         odampj = dabs(orbital_damping_factors(j))
      else
         odampj = 0._real64                           ! take the whole new orbital
      endif
      call dampor (j, inv, odampj)

!   Orthogonalize all orbitals of the same kappa in the order
!   fixed, spectroscopic, correlation orbitals. The order of
!   orbitals in the latter two classes are sorted according
!   to their self-consistency and energy.

      if (orthst) then
         !CALL orthor (J, inv)
         nwww = orbital_count
         call orthy (nwww, j, lsort)
      endif
!
!   Print details of iteration
!
      if (mpi_rank == 0)                                                &
         write (*, 302) orbital_principal_n(j),orbital_symmetry_labels(j),orbital_energy_parameters(j),method(j), &
            & orbital_origin_parameters(j),scnsty(j), &
!cjb DNORM-1 -> SQRT(DNORM)-1
!cjb                    DNORM - 1, ODAMPJ, JP, MF(J), INV, NNP
                        sqrt(dnorm)-1, odampj, jp, orbital_last_grid_index(j), inv, nnp
      dampmx = dmax1(dampmx,dabs(odampj))

  300 format(/,' Failure; equation for orbital ',1i2,1a2,&
         ' could not be solved using method ',1i1)
  301 format(/,/,' ****** Error in SUBROUTINE IMPROV ******'/,&
         ' Convergence not obtained'/)
  302 format (1x,1i2,1a2,1p,1d16.7,1x,1i2,d11.3,1d10.2,1d10.2,&
              0p,f6.3,1x,1i5,1x,1i5,1x,1i2,1x,1i2)
      return
      end subroutine improvmpi
