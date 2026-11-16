! MPI/serial migration of GRASP rdensity (upstream 488e875d).
! Accumulate one-body density matrices before contracting radial functions.
module density_calculation
    use vast_kind_param, only: double
    use parameter_def, only: nnnw
    use orb_c
    use wave_c
    use grid_c, only: n, r
    use eigv_c, only: evec
    use prnt_c, only: nvec, ivec
    use syma_c, only: iatjpo, iaspar
    use debug_c, only: cutoff
    use jlabl_c, only: jlbr, jlbp
    use density_parallel
    use natural_orbitals
    use onescalar_i
    use itjpo_i
    use ispar_i
    implicit none
    private
    public :: calculate_density
contains
    recursive subroutine calculate_density(name, nci, ordering)
        character(*), intent(in) :: name
        integer, intent(in) :: nci, ordering
        real(double), allocatable :: rho(:,:,:)
        real(double) :: angular(nnnw), coefficient
        integer, allocatable :: jvalue(:), parity(:)
        integer :: ic, ir, ia, ib, i, state, offset, different, delta

        allocate(rho(nw,nw,nvec), jvalue(ncf), parity(ncf))
        rho = 0.0d0
        do ic = 1,ncf
            jvalue(ic) = itjpo(ic)
            parity(ic) = ispar(ic)
        enddo
        ! Same cyclic column ownership as RHFS/RIS. Each process has its
        ! own angular-library workspace, including ranks with no columns.
        do ic = rank+1,ncf,workers
            if (rank == 0 .and. mod((ic-1)/workers,100) == 0) &
                write(*,'(a,i0)') 'Processing root CSF column: ', ic
            do ir = ic,ncf
                if (jvalue(ic) /= jvalue(ir) .or. parity(ic) /= parity(ir)) cycle
                ! Exact one-body selection rules, also used by RIS_MPI.
                different = 0
                do i = 1,nw
                    delta = int(iqa(i,ic))-int(iqa(i,ir))
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
                    do i = 1,nw
                        if (abs(angular(i)) <= cutoff) cycle
                        do state = 1,nvec
                            offset = (state-1)*ncf
                            coefficient = evec(ic+offset)*evec(ir+offset)*angular(i)
                            rho(i,i,state) = rho(i,i,state)+coefficient
                        enddo
                    enddo
                else
                    if (nak(ia) /= nak(ib) .or. abs(angular(1)) <= cutoff) cycle
                    do state = 1,nvec
                        offset = (state-1)*ncf
                        coefficient = evec(ic+offset)*evec(ir+offset)*angular(1)
                        rho(ia,ib,state) = rho(ia,ib,state)+coefficient
                        rho(ib,ia,state) = rho(ib,ia,state)+coefficient
                    enddo
                endif
            enddo
        enddo
        do state = 1,nvec
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
        real(double), intent(in) :: rho(:,:,:)
        real(double), allocatable :: radial(:)
        integer :: unit, ios, state, a, b, l
        character(256) :: filename

        filename = trim(name)//'.d'
        if (nci == 0) filename = trim(name)//'.cd'
        open(newunit=unit, file=filename, status='new', action='write', iostat=ios)
        if (ios /= 0) call density_fail('Cannot create '//trim(filename)//' (remove an existing output first)')
        allocate(radial(n))
        write(unit,*) '     r [au]          D(r)=r^2*rho(r)       rho(r)'
        write(unit,*)
        do state = 1,nvec
            radial = 0.0d0
            do a = 1,nw
                radial = radial+rho(a,a,state)*(pf(:n,a)**2+qf(:n,a)**2)
                do b = a+1,nw
                    if (nak(a) /= nak(b)) cycle
                    radial = radial+2.0d0*rho(a,b,state)*(pf(:n,a)*pf(:n,b)+qf(:n,a)*qf(:n,b))
                enddo
            enddo
            write(unit,'(1x,i3,5x,2a4)') ivec(state), jlbr(iatjpo(state)), jlbp((iaspar(state)+3)/2)
            do l = 2,n
                if (radial(l) > 0.0d0) write(unit,'(3es20.10e3)') r(l), radial(l), radial(l)/r(l)**2
            enddo
        enddo
        close(unit, iostat=ios)
        if (ios /= 0) call density_fail('Error writing '//filename)
    end subroutine

    recursive subroutine make_natural_orbitals(name, ordering, rho)
        character(*), intent(in) :: name
        integer, intent(in) :: ordering
        real(double), intent(in) :: rho(:,:,:)
        real(double), allocatable :: matrix(:,:), rotation(:,:), occupation(:)
        real(double), allocatable :: p(:,:), q(:,:), origin(:)
        integer, allocatable :: orbits(:)
        integer :: kappa, ndim, a, b, i, state, points, unit, ios, temp
        real(double) :: weight
        character(160) :: error

        weight = sum(real(iatjpo(:nvec), double))
        do kappa = minval(nak(:nw)),maxval(nak(:nw))
            ndim = count(nak(:nw) == kappa)
            if (ndim == 0) cycle
            allocate(orbits(ndim), matrix(ndim,ndim), rotation(ndim,ndim), occupation(ndim))
            orbits = pack([(i,i=1,nw)], nak(:nw) == kappa)
            ! Enumerate actual orbitals, so gaps in n and n>20 are supported.
            do a = 2,ndim
                b = a
                do while (b > 1)
                    if (np(orbits(b-1)) <= np(orbits(b))) exit
                    temp = orbits(b-1)
                    orbits(b-1) = orbits(b)
                    orbits(b) = temp
                    b = b-1
                enddo
            enddo
            matrix = 0.0d0
            do state = 1,nvec
                matrix = matrix+real(iatjpo(state),double)*rho(orbits,orbits,state)/weight
            enddo
            write(*,'(a,i5)') 'KAPPA: ', kappa
            write(*,'(a,i5)') 'NDIM: ', ndim
            write(*,*) 'density matrix:'
            do a = 1,ndim
                write(*,'(i5,*(1x,es24.16))') np(orbits(a)), matrix(:,a)
            enddo
            call diagonalize_density(matrix, ordering, occupation, rotation, error)
            if (len_trim(error) /= 0) call density_fail(trim(error))
            points = maxval(mf(orbits))
            allocate(p(points,ndim), q(points,ndim), origin(ndim))
            p = pf(:points,orbits)
            q = qf(:points,orbits)
            origin = pz(orbits)
            call transform_radial_block(rotation, p, q, origin)
            pf(:points,orbits) = p
            qf(:points,orbits) = q
            pz(orbits) = origin
            mf(orbits) = points
            write(*,*) 'occupations in output orbital order:'
            do a = 1,ndim
                write(*,'(i5,1x,es24.16)') np(orbits(a)), occupation(a)
            enddo
            write(*,*) 'eigenvectors in output orbital order (columns printed as rows):'
            do a = 1,ndim
                write(*,'(i5,*(1x,es24.16))') np(orbits(a)), rotation(:,a)
            enddo
            deallocate(orbits, matrix, rotation, occupation, p, q, origin)
        enddo
        open(newunit=unit, file=trim(name)//'.nw', form='unformatted', status='replace', iostat=ios)
        if (ios /= 0) call density_fail('Cannot create natural-orbital file')
        write(unit) 'G92RWF'
        do i = 1,nw
            ! E is retained for GRASP file compatibility; it is not an NO energy.
            write(unit) np(i), nak(i), e(i), mf(i)
            write(unit) pz(i), pf(:mf(i),i), qf(:mf(i),i)
            write(unit) r(:mf(i))
        enddo
        close(unit, iostat=ios)
        if (ios /= 0) call density_fail('Error writing natural-orbital file')
        write(*,*) 'Natural orbitals written to ', trim(name)//'.nw'
    end subroutine
end module density_calculation
