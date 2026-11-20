! Compile the unmodified production MATRIX body; only assembly/eigensolution
! and memory-record input are controlled seams. BLAS is the linked real library.
module ci_fixture
    use eigenstate_expansions, only: evec => csf_mixing_coefficients, eval => shifted_eigenvalues
    implicit none (type, external)
    integer :: solves = 0
contains
    subroutine eigenpairs()
        use, intrinsic :: iso_fortran_env, only: real64
        solves = solves + 1
        eval = [-2.0_real64, -1.0_real64]
        if (solves == 1) then
            evec(3:6) = [1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64]
        else
            evec(3:6) = [0.6_real64, 0.8_real64, -0.8_real64, 0.6_real64]
        end if
end subroutine eigenpairs
end module ci_fixture

program ci_damping
    use ci_fixture, only: eigenpairs, eval, evec, solves
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use scf_damping_factors, only: cdamp => mixing_damping_factors
    use atomic_calculation_state, only: iccmin, ncmin, ncmax
    use debug_output_settings, only: ldbpg
    use hamiltonian_symmetry_blocks, only: nblock => symmetry_block_count, ncfblk => csf_counts_by_block, nevblk => &
       & eigenvector_counts_by_block, ncmaxblk
    use sparse_hamiltonian_storage, only: &
        emt => hamiltonian_values, irow => matrix_partner_indices, iendc => stored_segment_end_offsets
    use mcp_file_metadata, only: kmaxf => max_mcp_file_multipole_rank, mcp_input_units
    use calculation_file_units, only: mixing_output_unit
    use orbital_csf_basis, only: ncf => csf_count, nw => orbital_count
    use symmetry_block_offsets, only: ncftot => total_csf_count, nvecsiz => mixing_coefficient_storage_size, ncfpast => &
       & csf_offsets_by_block, ncminpast => eigenvalue_offsets_by_block, nevecpast => mixing_offsets_by_block
    use block_mean_energies, only: eavblk
    use level_symmetry_metadata, only: jpgg
#ifdef CI_MPI
    use mpi_state, only: host => processor_name, ierr => mpi_error_code, istat => mpi_status_buffer, &
       & lenhost => processor_name_length, myid => mpi_rank, nprocs => mpi_process_count, mpi_finalize, &
       & startmpi
    use matrixmpi_I, only: matrixmpi
#else
    use mpi_s, only: myid => mpi_rank, nprocs => mpi_process_count
    use matrix_I, only: matrix
#endif
#ifdef CI_MEM
    use rmcdhf_mem_C, only: imin1_30, imax1_30, imin2_30, imax2_30
#endif
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)
    integer :: channel_suffix, i, records, local_size
    real(real64), parameter :: tol = 1e-12_real64
#ifdef CI_MPI
    call startmpi(myid,nprocs,host,lenhost)
#else
    myid = 0
    nprocs = 1
#endif
    nblock=1; ncf=2; ncftot=2; nw=1; ncmin=2; ncmax=2; nvecsiz=8; kmaxf=1
    ldbpg=.false.; jpgg=1
    allocate(ncfblk(1), nevblk(1), ncmaxblk(1), ncfpast(1), ncminpast(1), nevecpast(1), eavblk(1))
    allocate(iccmin(2), cdamp(2), evec(8), eval(2))
    ncfblk=2; nevblk=2; ncmaxblk=2; ncfpast=0; ncminpast=0; nevecpast=2
    iccmin=[1,2]; cdamp=[0.0_real64,0.5_real64]; evec=99.0_real64
    nullify(emt,irow,iendc)
    local_size = count([(mod(i-1,nprocs)==myid, i=1,2)])
#ifdef CI_MEM
    allocate(imin1_30(1),imax1_30(1),imin2_30(1),imax2_30(1))
    imin1_30=0; imax1_30=2; imin2_30=0; imax2_30=local_size
#else
    allocate(mcp_input_units(30:33))
    do channel_suffix=30,33
        open(newunit=mcp_input_units(channel_suffix),status="scratch",form="unformatted")
        records=3
        if(channel_suffix==30) records=6
        do i=1,records
            write(mcp_input_units(channel_suffix)) 0
        end do
        if(channel_suffix==30) then
            write(mcp_input_units(channel_suffix)) "MCP",1,2
            write(mcp_input_units(channel_suffix)) local_size
            write(mcp_input_units(channel_suffix)) [(i,i=1,local_size)], [(i,i=1,local_size)]
        end if
    end do
#endif
    if(myid==0) then
        open(newunit=mixing_output_unit,status="scratch",form="unformatted")
        write(mixing_output_unit) "G92MIX"
        write(mixing_output_unit) 2,2,1,2,8,1
    end if
    call matrix_step(.true.)
    ! With fcheck=all, the production allocation/zeroing must include IENDC(0).
    if(maxval(abs(evec(3:6)-[1.0_real64,0.0_real64,0.0_real64,1.0_real64]))>tol) &
        error stop "FIRST must bypass damping"
    call matrix_step(.false.)
    if(.not.all(ieee_is_finite(evec))) error stop "nonfinite damped roots"
    if(maxval(abs(evec(3:6)-[0.6_real64,0.8_real64,0.8_real64,-0.6_real64]))>tol) &
        error stop "DDOT damping branch differs from analytic roots"
    if(abs(dot_product(evec(3:4),evec(5:6)))>tol) error stop "orthogonality"
    if(abs(sum(evec(3:4)**2)-1.0_real64)>tol.or.abs(sum(evec(5:6)**2)-1.0_real64)>tol) &
        error stop "normalization"
    if(any(evec(1:2)/=99.0_real64).or.any(evec(7:8)/=99.0_real64)) error stop "block offset guards"
    cdamp=0.0_real64
    call matrix_step(.false.)
    if(maxval(abs(evec(3:6)-[0.6_real64,0.8_real64,-0.8_real64,0.6_real64]))>tol) &
        error stop "zero damping control"
    if(solves/=3) error stop "fixture did not execute three production calls"
    if(myid==0) then
        close(mixing_output_unit)
        print *, "PASS: production FIRST / nonzero multi-root DDOT / zero damping"
    end if
#ifndef CI_MEM
    do channel_suffix=30,33
        close(mcp_input_units(channel_suffix))
    end do
#endif
#ifdef CI_MPI
    call MPI_Finalize(ierr)
#endif
contains
    subroutine matrix_step(first)
        logical,intent(in)::first
#ifdef CI_MPI
        call matrixmpi(first)
#else
        call matrix(first)
#endif
end subroutine matrix_step
end program ci_damping

#ifdef CI_MPI
subroutine maneigmpi(dvdfirst,lprint,jblock,ncfpat,ncminpat,nevecpat,ncftot)
#else
subroutine maneig(dvdfirst,lprint,jblock,ncfpat,ncminpat,nevecpat,ncftot)
#endif
    use ci_fixture
    implicit none (type, external)
    logical,intent(in)::dvdfirst
    logical::lprint
    integer,intent(in)::jblock,ncfpat,ncminpat,nevecpat
    integer::ncftot
    if(jblock/=1.or.nevecpat/=2.or.ncftot/=2) error stop "eigensolver seam arguments"
    call eigenpairs()
#ifdef CI_MPI
end subroutine maneigmpi
#else
end subroutine maneig
#endif

#ifdef CI_MEM
subroutine setham(jblock,myid,nprocs,ncontr_tot,ncontr_tot2)
#else
subroutine setham(jblock,myid,nprocs)
#endif
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use sparse_hamiltonian_storage, only: emt => hamiltonian_values
    implicit none (type, external)
    integer,intent(in)::jblock,myid,nprocs
#ifdef CI_MEM
    integer(int64)::ncontr_tot(20),ncontr_tot2(20)
#endif
    if(jblock/=1) error stop "Hamiltonian seam block"
    emt=0.0_real64
end subroutine setham

#ifdef CI_MEM
subroutine read03_mem(jblock)
    use sparse_hamiltonian_storage, only: irow => matrix_partner_indices, iendc => stored_segment_end_offsets
#ifdef CI_MPI
    use mpi_state, only: myid => mpi_rank, nprocs => mpi_process_count
#else
    use mpi_s, only: myid => mpi_rank, nprocs => mpi_process_count
#endif
    implicit none (type, external)
    integer,intent(in)::jblock
    integer::i,k
    if(jblock/=1) error stop "memory seam block"
    k=0
    do i=myid+1,2,nprocs
        k=k+1
        iendc(i)=k
        irow(k)=i
    end do
end subroutine read03_mem
subroutine read1_mem(nfile,ncontr_tot2,lab,ncontr)
    use, intrinsic :: iso_fortran_env, only: int64
    implicit none (type, external)
    integer,intent(in)::nfile
    integer(int64),intent(in)::ncontr_tot2
    integer,intent(out)::lab,ncontr
    error stop "unexpected empty block seam"
end subroutine read1_mem
#endif
