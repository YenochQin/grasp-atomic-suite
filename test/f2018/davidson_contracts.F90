! Production initial estimates, sparse multiplication and Davidson solver callback contract.
program davidson_contracts
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use real_vector_initialization, only: initialize_real_vector
    use davidson_sparse_initial_estimates, only: estimate_davidson_sparse_eigenpairs
    use sparse_eigensolver_initial_basis, only: estimate_sparse_initial_basis
    use davidson_eigensolver, only: solve_davidson_eigenpairs
    use sparse_symmetric_matvec, only: multiply_sparse_symmetric_matrix
    use davidson_projected_matrix, only: update_davidson_projected_matrix
    use davidson_residual_expansion, only: expand_davidson_residuals
    use sparse_hamiltonian_storage, only: hamiltonian_values, matrix_partner_indices, stored_segment_end_offsets
    implicit none (type, external)
    integer, parameter :: n=3, nd=6, lim=6, lwork=2*nd*lim+2*lim*lim+11*lim+1, liwork=6*lim+1
    integer :: niv, ilow, ihigh, mblock, maxiter, loops, nmv, ierr, i
    integer :: jcol(0:3), rows(3), blocks(3), selected(lim), iwork(liwork)
    real(real64) :: values(3), basis(8), work(lwork), diag(nd), h(n,n), big(nd,nd)
    real(real64) :: result(n,2), vectors(n,2), crite, critc, critr, guard(7)
    logical :: hiend
    character(16) :: scenario
    real(real64) :: expected
    guard = -99
    call initialize_real_vector(3, 2.0_real64, guard(2:6), 2)
    if (any(guard /= [-99.0_real64,2.0_real64,-99.0_real64,2.0_real64, &
                     -99.0_real64,2.0_real64,-99.0_real64])) error stop "DINIT stride/guards"
    call initialize_real_vector(3, 4.0_real64, guard(2:6), -2)
    if (any(guard(2:6:2) /= 4.0_real64)) error stop "DINIT negative stride"
    call initialize_real_vector(0, 7.0_real64, guard, 1)
    if (guard(1) /= -99.0_real64) error stop "DINIT zero count"

    values = [1.0_real64,2.0_real64,4.0_real64]
    rows = [1,2,3]
    jcol = [0,1,2,3]
    blocks = 1
    niv = 2
    basis = -99
    call estimate_davidson_sparse_eigenpairs(n,n,niv,values,jcol,rows,basis,blocks,1)
    call check_estimates()
    basis = -99
    call estimate_sparse_initial_basis(n,n,niv,basis,values,int(jcol,int64),rows)
    call check_estimates()

    h = reshape([2.0_real64,1.0_real64,0.0_real64, &
                 1.0_real64,2.0_real64,0.0_real64, &
                 0.0_real64,0.0_real64,4.0_real64], [n,n])
    allocate(hamiltonian_values(4),matrix_partner_indices(4),stored_segment_end_offsets(0:n))
    hamiltonian_values = [2.0_real64,1.0_real64,2.0_real64,4.0_real64]
    matrix_partner_indices = [1,1,2,3]
    stored_segment_end_offsets = [0,1,3,4]
    vectors = reshape([1.0_real64,2.0_real64,3.0_real64,4.0_real64,5.0_real64,6.0_real64],[n,2])
    result = -99
    call multiply_sparse_symmetric_matrix(n,2,vectors,result)
    if (maxval(abs(result-matmul(h,vectors))) > 1e-12_real64) error stop "SPICMV2 sparse result"
    deallocate(hamiltonian_values,matrix_partner_indices,stored_segment_end_offsets)
    call check_partial_updates()

    ! Davidson solver reserves three restart slots, so use LIM=6 even for this small block.
    call get_command_argument(1, scenario)
    if (scenario == "cluster") then
        call check_cluster()
        stop
    end if
    big=0
    select case (trim(scenario))
    case ("")
        do i=1,nd
            big(i,i)=2.0_real64
        end do
        do i=1,nd-1
            big(i,i+1)=-1.0_real64
            big(i+1,i)=-1.0_real64
        end do
        expected=2.0_real64-2.0_real64*cos(acos(-1.0_real64)/real(nd+1,real64))
    case ("early")
        ! An invariant two-dimensional block reaches an exact Ritz pair
        ! before a third expansion. It must never normalize a zero residual.
        do i=1,nd
            big(i,i)=real(i+1,real64)
        end do
        big(2,2)=2.0_real64
        big(1,2)=1.0_real64
        big(2,1)=1.0_real64
        expected=1.0_real64
    case default
        error stop "Unknown Davidson scenario"
    end select
    work=0
    work(1)=1
    work(nd+1)=big(1,1)
    do i=1,nd
        diag(i)=big(i,i)
    end do
    niv=1
    ilow=1
    ihigh=1
    mblock=1
    maxiter=30
    selected=0
    iwork=0
    loops=0
    hiend=.false.
    crite=1e-14_real64
    critc=1e-12_real64
    critr=1e-12_real64
    call solve_davidson_eigenpairs(multiply,nd,lim,diag,ilow,ihigh,selected,niv,mblock,crite,critc,critr, &
              maxiter,work,lwork,iwork,liwork,hiend,loops,nmv,ierr)
    if (ierr /= 0 .or. nmv < 1) error stop "Davidson solver callback execution"
    if (.not. all(ieee_is_finite(work(:nd+1)))) error stop "Davidson solver nonfinite eigenpair"
    if (abs(work(nd+1)-expected) > 1e-12_real64) &
        error stop "Davidson solver eigenvalue"
    if (maxval(abs(matmul(big,work(:nd))-work(nd+1)*work(:nd))) > 1e-12_real64) error stop "Davidson solver residual"
    if (abs(sum(work(:nd)**2)-1.0_real64) > 1e-12_real64) error stop "Davidson solver normalization"
    print *, "PASS Davidson interfaces and numerical contracts"
contains
    subroutine check_cluster()
        integer, parameter :: chain = 7, order = 2*chain, limit = 8, roots = 2
        integer, parameter :: real_size = 2*order*limit+2*limit*limit+11*limit+roots
        integer, parameter :: integer_size = 6*limit+roots
        real(real64), parameter :: gap = 1.0e-7_real64, tolerance = 1.0e-10_real64
        real(real64) :: diagonal(order), matrix(order,order), reference(order,roots), actual(order,roots)
        real(real64) :: workspace(real_size), expected_values(roots), actual_values(roots), gram(roots,roots)
        real(real64) :: projector(order,order), residual(order,roots), angle
        integer :: index, root, column, partner, position, iteration, count, block_size, attempt, iteration_limit
        integer :: initial_vector_count
        integer :: first_root, last_root, status, selection(limit), scratch(integer_size)
        logical :: high_end
        real(real64) :: energy_tol, coefficient_tol, residual_tol

        ! Two identical tridiagonal chains, with a known small shift on the second chain.
        ! A trial space narrower than the matrix order forces a restart in the cold solve.
        angle = acos(-1.0_real64)/(chain+1)
        expected_values(1) = 2.0_real64 - 2.0_real64*cos(angle)
        expected_values(2) = expected_values(1) + gap
        matrix = 0.0_real64
        reference = 0.0_real64
        position = 0
        allocate(hamiltonian_values(2*order-2),matrix_partner_indices(2*order-2),stored_segment_end_offsets(0:order))
        stored_segment_end_offsets(0) = 0_int64
        do column = 1, order
            index = mod(column-1,chain)+1
            root = 1+(column-1)/chain
            diagonal(column) = 2.0_real64 + (root-1)*gap
            matrix(column,column) = diagonal(column)
            reference(column,root) = sqrt(2.0_real64/(chain+1))*sin(index*angle)
            if (index > 1) then
                partner = column-1
                matrix(partner,column) = -1.0_real64
                matrix(column,partner) = -1.0_real64
                position = position+1
                hamiltonian_values(position) = -1.0_real64
                matrix_partner_indices(position) = partner
            end if
            position = position+1
            hamiltonian_values(position) = diagonal(column)
            matrix_partner_indices(position) = column
            stored_segment_end_offsets(column) = position
        end do
        projector = matmul(reference,transpose(reference))
        workspace = 0.0_real64
        workspace(1) = 1.0_real64
        workspace(order+chain+1) = 1.0_real64

        do attempt = 1, 2
            ! The second solve reuses the two returned vectors with fresh scratch space.
            workspace(order*roots+1:) = 0.0_real64
            scratch = 0
            selection = 0
            first_root = 1
            last_root = roots
            block_size = 2
            initial_vector_count = roots
            if (attempt == 2) then
                initial_vector_count = roots+1
                workspace(order*roots+1:order*roots+chain) = &
                    sqrt(2.0_real64/(chain+1))*sin(2*[(index,index=1,chain)]*angle)
            end if
            iteration_limit = 500
            iteration = 0
            high_end = .false.
            energy_tol = 1.0e-14_real64
            coefficient_tol = 1.0e-12_real64
            residual_tol = 1.0e-12_real64
            call solve_davidson_eigenpairs(multiply_sparse_symmetric_matrix,order,limit,diagonal,first_root,last_root, &
                selection,initial_vector_count,block_size,energy_tol,coefficient_tol,residual_tol,iteration_limit, &
                workspace,real_size,scratch,integer_size,high_end,iteration,count,status)
            if (status /= 0 .or. count < 1) error stop "Davidson cluster solve failed"
            if (attempt == 1 .and. count <= limit) error stop "cold cluster did not exercise a restart"
            if (attempt == 2 .and. count /= roots+1) error stop "cluster warm start expanded an already converged basis"
            if (.not. all(ieee_is_finite(workspace(:order*roots+roots)))) error stop "nonfinite cluster result"
            actual_values = workspace(order*roots+1:order*roots+roots)
            actual = reshape(workspace(:order*roots),shape(actual))
            if (maxval(abs(actual_values-expected_values)) > tolerance) error stop "wrong cluster eigenvalues"
            residual = matmul(matrix,actual) - actual*spread(actual_values,1,order)
            write (*,"(a,i0,a,i0,a,2es15.6)") "Cluster attempt ", attempt, ", products ", count, &
                ", residuals ", sqrt(sum(residual**2,dim=1))
            if (maxval(sqrt(sum(residual**2,dim=1))) > tolerance) error stop "cluster residual"
            gram = matmul(transpose(actual),actual)
            do root = 1, roots
                gram(root,root) = gram(root,root)-1.0_real64
            end do
            if (maxval(abs(gram)) > tolerance) error stop "cluster orthogonality"
            if (maxval(abs(matmul(actual,transpose(actual))-projector)) > tolerance) error stop "cluster subspace"
        end do
        deallocate(hamiltonian_values,matrix_partner_indices,stored_segment_end_offsets)
        print *, "PASS Davidson near-degenerate roots, restart and reusable vectors"
    end subroutine check_cluster

    subroutine check_partial_updates()
        real(real64) :: b(3,3), ab(3,3), packed(6), coeff(6), eigen(3), old(2), scratch(3), shifts(3)
        integer :: flags(2), pending(2), count_new
        logical :: done
        b=0
        b(1,1)=1; b(2,2)=1; b(3,3)=1
        packed=-99
        call update_davidson_projected_matrix(3,3,0,1,b,h,packed)
        if(abs(packed(1)-2.0_real64)>0.or.maxval(abs(packed(2:)+99.0_real64))>0) &
            error stop "ADDS first-column guards"
        call update_davidson_projected_matrix(3,3,1,1,b,h,packed)
        call update_davidson_projected_matrix(3,3,2,1,b,h,packed)
        if(.not.all(ieee_is_finite(packed))) error stop "Nonfinite projected matrix"
        if(maxval(abs(packed-[2.0_real64,1.0_real64,2.0_real64,0.0_real64,0.0_real64,4.0_real64]))>0) &
            error stop "ADDS must preserve previous projected columns"
        ! Root 2 is already converged. Updating root 1 must retain that flag.
        ab=0; ab(1,1)=1; ab(2,2)=2
        coeff=[1.0_real64,0.0_real64,0.0_real64,1.0_real64,0.0_real64,0.0_real64]
        eigen=[1.0_real64,2.0_real64,3.0_real64];old=[0.0_real64,0.0_real64]
        flags=[0,1];pending=[1,2];count_new=1;scratch=-99;shifts=-99
        call expand_davidson_residuals(3,2,3,1,2,1e-12_real64,count_new,pending,coeff,eigen,old,ab,b,flags,scratch,shifts,done)
        if(.not.done.or.any(flags/=1)) error stop "NEWVEC must retain other roots convergence state"
    end subroutine check_partial_updates
    subroutine check_estimates()
        if (any(abs(basis(7:8)-values(:2)) > 1e-12_real64)) error stop "INIEST eigenvalues"
        if (abs(abs(basis(1))-1.0_real64) > 1e-12_real64) error stop "INIEST first vector"
        if (abs(abs(basis(5))-1.0_real64) > 1e-12_real64) error stop "INIEST second vector"
        if (any(abs(basis([2,3,4,6])) > 1e-12_real64)) error stop "INIEST zero components"
    end subroutine check_estimates
    subroutine multiply(size_n, count, b, c)
        integer, intent(in) :: size_n, count
        real(real64), intent(in) :: b(size_n,count)
        real(real64), intent(inout) :: c(size_n,count)
        c = matmul(big,b)
    end subroutine multiply
end program davidson_contracts
