! MPI/serial migration of GRASP rdensity (upstream 488e875d).
! Accumulate one-body density matrices before contracting radial functions.
module density_calculation
    use vast_kind_param, only: DOUBLE
    use parameter_def, only: NNNW
    use orb_C
    use wave_C
    use grid_C, only: N, R
    use eigv_C, only: EVEC
    use prnt_C, only: NVEC, IVEC
    use syma_C, only: IATJPO, IASPAR
    use debug_C, only: CUTOFF
    use jlabl_C, only: JLBR, JLBP
    use density_parallel
    use natural_orbitals
    use onescalar_I
    use itjpo_I
    use ispar_I
    implicit none
    private
    public :: calculate_density
contains
    recursive subroutine calculate_density(name, nci, ordering)
        character(*), intent(in) :: name
        integer, intent(in) :: nci, ordering
        real(DOUBLE), allocatable :: rho(:,:,:)
        real(DOUBLE) :: angular(NNNW), coefficient
        integer, allocatable :: jvalue(:), parity(:)
        integer :: ic, ir, ia, ib, i, state, offset, different, delta

        allocate(rho(NW,NW,NVEC), jvalue(NCF), parity(NCF))
        rho = 0.0D0
        do ic = 1,NCF
            jvalue(ic) = itjpo(ic)
            parity(ic) = ispar(ic)
        enddo
        ! Same cyclic column ownership as RHFS/RIS. Each process has its
        ! own angular-library workspace, including ranks with no columns.
        do ic = rank+1,NCF,workers
            if (rank == 0 .and. mod((ic-1)/workers,100) == 0) &
                write(*,'(a,i0)') 'Processing root CSF column: ', ic
            do ir = ic,NCF
                if (jvalue(ic) /= jvalue(ir) .or. parity(ic) /= parity(ir)) cycle
                ! Exact one-body selection rules, also used by RIS_MPI.
                different = 0
                do i = 1,NW
                    delta = int(IQA(i,ic))-int(IQA(i,ir))
                    if (delta == 0) cycle
                    different = different+1
                    if (abs(delta) > 1 .or. different > 2) then
                        different = 3
                        exit
                    endif
                enddo
                if (different /= 0 .and. different /= 2) cycle
                if (different == 0 .and. ic /= ir) cycle
                call onescalar(ic, ir, ia, ib, angular)
                if (ia == 0) cycle
                if (ia == ib) then
                    do i = 1,NW
                        if (abs(angular(i)) <= CUTOFF) cycle
                        do state = 1,NVEC
                            offset = (state-1)*NCF
                            coefficient = EVEC(ic+offset)*EVEC(ir+offset)*angular(i)
                            rho(i,i,state) = rho(i,i,state)+coefficient
                        enddo
                    enddo
                else
                    if (NAK(ia) /= NAK(ib) .or. abs(angular(1)) <= CUTOFF) cycle
                    do state = 1,NVEC
                        offset = (state-1)*NCF
                        coefficient = EVEC(ic+offset)*EVEC(ir+offset)*angular(1)
                        rho(ia,ib,state) = rho(ia,ib,state)+coefficient
                        rho(ib,ia,state) = rho(ib,ia,state)+coefficient
                    enddo
                endif
            enddo
        enddo
        do state = 1,NVEC
            call sum_density(rho(:,:,state))
        enddo
        if (rank == 0) then
            ! Use the original radial basis before rotating PF and QF.
            call write_density(name, nci, rho)
            call make_natural_orbitals(name, ordering, rho)
        endif
    end subroutine

    recursive subroutine write_density(name, nci, rho)
        character(*), intent(in) :: name
        integer, intent(in) :: nci
        real(DOUBLE), intent(in) :: rho(:,:,:)
        real(DOUBLE), allocatable :: radial(:)
        integer :: unit, ios, state, a, b, l
        character(256) :: filename

        filename = trim(name)//'.d'
        if (nci == 0) filename = trim(name)//'.cd'
        open(newunit=unit, file=filename, status='new', action='write', iostat=ios)
        if (ios /= 0) call density_fail('Cannot create '//trim(filename)//' (remove an existing output first)')
        allocate(radial(N))
        write(unit,*) '     r [au]          D(r)=r^2*rho(r)       rho(r)'
        write(unit,*)
        do state = 1,NVEC
            radial = 0.0D0
            do a = 1,NW
                radial = radial+rho(a,a,state)*(PF(:N,a)**2+QF(:N,a)**2)
                do b = a+1,NW
                    if (NAK(a) /= NAK(b)) cycle
                    radial = radial+2.0D0*rho(a,b,state)*(PF(:N,a)*PF(:N,b)+QF(:N,a)*QF(:N,b))
                enddo
            enddo
            write(unit,'(1x,i3,5x,2a4)') IVEC(state), JLBR(IATJPO(state)), JLBP((IASPAR(state)+3)/2)
            do l = 2,N
                if (radial(l) > 0.0D0) write(unit,'(3es20.10e3)') R(l), radial(l), radial(l)/R(l)**2
            enddo
        enddo
        close(unit, iostat=ios)
        if (ios /= 0) call density_fail('Error writing '//filename)
    end subroutine

    recursive subroutine make_natural_orbitals(name, ordering, rho)
        character(*), intent(in) :: name
        integer, intent(in) :: ordering
        real(DOUBLE), intent(in) :: rho(:,:,:)
        real(DOUBLE), allocatable :: matrix(:,:), rotation(:,:), occupation(:)
        real(DOUBLE), allocatable :: p(:,:), q(:,:), origin(:)
        integer, allocatable :: orbits(:)
        integer :: kappa, ndim, a, b, i, state, points, unit, ios, temp
        real(DOUBLE) :: weight
        character(160) :: error

        weight = sum(real(IATJPO(:NVEC), DOUBLE))
        do kappa = minval(NAK(:NW)),maxval(NAK(:NW))
            ndim = count(NAK(:NW) == kappa)
            if (ndim == 0) cycle
            allocate(orbits(ndim), matrix(ndim,ndim), rotation(ndim,ndim), occupation(ndim))
            orbits = pack([(i,i=1,NW)], NAK(:NW) == kappa)
            ! Enumerate actual orbitals, so gaps in n and n>20 are supported.
            do a = 2,ndim
                b = a
                do while (b > 1)
                    if (NP(orbits(b-1)) <= NP(orbits(b))) exit
                    temp = orbits(b-1)
                    orbits(b-1) = orbits(b)
                    orbits(b) = temp
                    b = b-1
                enddo
            enddo
            matrix = 0.0D0
            do state = 1,NVEC
                matrix = matrix+real(IATJPO(state),DOUBLE)*rho(orbits,orbits,state)/weight
            enddo
            write(*,'(a,i5)') 'KAPPA: ', kappa
            write(*,'(a,i5)') 'NDIM: ', ndim
            write(*,*) 'density matrix:'
            do a = 1,ndim
                write(*,'(i5,*(1x,es24.16))') NP(orbits(a)), matrix(:,a)
            enddo
            call diagonalize_density(matrix, ordering, occupation, rotation, error)
            if (len_trim(error) /= 0) call density_fail(trim(error))
            points = maxval(MF(orbits))
            allocate(p(points,ndim), q(points,ndim), origin(ndim))
            p = PF(:points,orbits)
            q = QF(:points,orbits)
            origin = PZ(orbits)
            call transform_radial_block(rotation, p, q, origin)
            PF(:points,orbits) = p
            QF(:points,orbits) = q
            PZ(orbits) = origin
            MF(orbits) = points
            write(*,*) 'occupations in output orbital order:'
            do a = 1,ndim
                write(*,'(i5,1x,es24.16)') NP(orbits(a)), occupation(a)
            enddo
            write(*,*) 'eigenvectors in output orbital order (columns printed as rows):'
            do a = 1,ndim
                write(*,'(i5,*(1x,es24.16))') NP(orbits(a)), rotation(:,a)
            enddo
            deallocate(orbits, matrix, rotation, occupation, p, q, origin)
        enddo
        open(newunit=unit, file=trim(name)//'.nw', form='unformatted', status='replace', iostat=ios)
        if (ios /= 0) call density_fail('Cannot create natural-orbital file')
        write(unit) 'G92RWF'
        do i = 1,NW
            ! E is retained for GRASP file compatibility; it is not an NO energy.
            write(unit) NP(i), NAK(i), E(i), MF(i)
            write(unit) PZ(i), PF(:MF(i),i), QF(:MF(i),i)
            write(unit) R(:MF(i))
        enddo
        close(unit, iostat=ios)
        if (ios /= 0) call density_fail('Error writing natural-orbital file')
        write(*,*) 'Natural orbitals written to ', trim(name)//'.nw'
    end subroutine
end module density_calculation
