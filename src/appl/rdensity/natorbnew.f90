! MPI/serial migration of GRASP rdensity (upstream 488e875d).
! Accumulate one-body density matrices before contracting radial functions.
module density_calculation
    use parameter_def, only: max_orbitals
    use orb_c, only: csf_count, csf_subshell_occupations, double, max_orbitals, orbital_count, &
        orbital_energy_parameters, orbital_kappa, orbital_principal_n
    use wave_c, only: double, max_orbitals, orbital_last_grid_index, orbital_origin_parameters, mfff, pzff, mfii, &
        pzii, radial_large_components, radial_small_components, pfff, qfff, pfii, qfii
    use grid_c, only: radial_point_count, radial_grid
    use eigv_c, only: csf_mixing_coefficients
    use prnt_c, only: nvec, ivec
    use syma_c, only: iatjpo, iaspar
    use debug_c, only: cutoff
    use jlabl_c, only: jlbr, jlbp
    use density_parallel, only: density_fail, rank, sum_density, workers
    use natural_orbitals, only: diagonalize_density, transform_radial_block
    use onescalar_i, only: onescalar
    use itjpo_i, only: itjpo
    use ispar_i, only: ispar
    implicit none (type, external)
    private
    public :: calculate_density
contains
    recursive subroutine calculate_density(name, nci, ordering)
        use, intrinsic :: iso_fortran_env, only: real64
        character(*), intent(in) :: name
        integer, intent(in) :: nci, ordering
        real(real64), allocatable :: rho(:,:,:)
        real(real64) :: angular(max_orbitals), coefficient
        integer, allocatable :: jvalue(:), parity(:)
        integer :: ic, ir, ia, ib, i, state, offset, different, delta

        allocate(rho(orbital_count,orbital_count,nvec), jvalue(csf_count), parity(csf_count))
        rho = 0.0_real64
        do ic = 1,csf_count
            jvalue(ic) = itjpo(ic)
            parity(ic) = ispar(ic)
        end do
        ! Same cyclic column ownership as RHFS/RIS. Each process has its
        ! own angular-library workspace, including ranks with no columns.
        do ic = rank+1,csf_count,workers
            if (rank == 0 .and. mod((ic-1)/workers,100) == 0) &
                write(*,"(a,i0)") "Processing root CSF column: ", ic
            do ir = ic,csf_count
                if (jvalue(ic) /= jvalue(ir) .or. parity(ic) /= parity(ir)) cycle
                ! Exact one-body selection rules, also used by RIS_MPI.
                different = 0
                do i = 1,orbital_count
                    delta = int(csf_subshell_occupations(i,ic))-int(csf_subshell_occupations(i,ir))
                    if (delta == 0) cycle
                    different = different+1
                    if (abs(delta) > 1 .or. different > 2) then
                        different = 3
                        exit
                    end if
                end do
                if (different /= 0 .and. different /= 2) cycle
                if (different == 0 .and. ic /= ir) cycle
                call onescalar(ic, ir, ia, ib, angular)
                if (ia == 0) cycle
                if (ia == ib) then
                    do i = 1,orbital_count
                        if (abs(angular(i)) <= cutoff) cycle
                        do state = 1,nvec
                            offset = (state-1)*csf_count
                            coefficient = csf_mixing_coefficients(ic+offset)*csf_mixing_coefficients(ir+offset)*angular(i)
                            rho(i,i,state) = rho(i,i,state)+coefficient
                        end do
                    end do
                else
                    if (orbital_kappa(ia) /= orbital_kappa(ib) .or. abs(angular(1)) <= cutoff) cycle
                    do state = 1,nvec
                        offset = (state-1)*csf_count
                        coefficient = csf_mixing_coefficients(ic+offset)*csf_mixing_coefficients(ir+offset)*angular(1)
                        rho(ia,ib,state) = rho(ia,ib,state)+coefficient
                        rho(ib,ia,state) = rho(ib,ia,state)+coefficient
                    end do
                end if
            end do
        end do
        do state = 1,nvec
            call sum_density(rho(:,:,state))
        end do
        if (rank == 0) then
            ! Use the original radial basis before rotating PF and QF.
            call write_density(name, nci, rho)
            call make_natural_orbitals(name, ordering, rho)
        end if
    end subroutine calculate_density

    recursive subroutine write_density(name, nci, rho)
        use, intrinsic :: iso_fortran_env, only: real64
        character(*), intent(in) :: name
        integer, intent(in) :: nci
        real(real64), intent(in) :: rho(:,:,:)
        real(real64), allocatable :: radial(:)
        integer :: unit, ios, state, a, b, l
        character(256) :: filename

        filename = trim(name)//".d"
        if (nci == 0) filename = trim(name)//".cd"
        open(newunit=unit, file=filename, status="new", action="write", iostat=ios)
        if (ios /= 0) call density_fail("Cannot create "//trim(filename)//" (remove an existing output first)")
        allocate(radial(radial_point_count))
        write(unit,*) "     r [au]          D(r)=r^2*rho(r)       rho(r)"
        write(unit,*)
        do state = 1,nvec
            radial = 0.0_real64
            do a = 1,orbital_count
                radial = radial+rho(a,a,state)*(radial_large_components(:radial_point_count,a)**2+ &
                   & radial_small_components(:radial_point_count,a)**2)
                do b = a+1,orbital_count
                    if (orbital_kappa(a) /= orbital_kappa(b)) cycle
                    radial = radial+2.0_real64*rho(a,b,state)*(radial_large_components(:radial_point_count,a)* &
                       & radial_large_components(:radial_point_count,b)+radial_small_components(:radial_point_count,a)* &
                       & radial_small_components(:radial_point_count,b))
                end do
            end do
            write(unit,"(1x,i3,5x,2a4)") ivec(state), jlbr(iatjpo(state)), jlbp((iaspar(state)+3)/2)
            do l = 2,radial_point_count
                if (radial(l) > 0.0_real64) write(unit,"(3es20.10e3)") radial_grid(l), radial(l), radial(l)/radial_grid(l)**2
            end do
        end do
        close(unit, iostat=ios)
        if (ios /= 0) call density_fail("Error writing "//filename)
    end subroutine write_density

    recursive subroutine make_natural_orbitals(name, ordering, rho)
        use, intrinsic :: iso_fortran_env, only: real64
        character(*), intent(in) :: name
        integer, intent(in) :: ordering
        real(real64), intent(in) :: rho(:,:,:)
        real(real64), allocatable :: matrix(:,:), rotation(:,:), occupation(:)
        real(real64), allocatable :: p(:,:), q(:,:), origin(:)
        integer, allocatable :: orbits(:)
        integer :: kappa, ndim, a, b, i, state, points, unit, ios, temp
        real(real64) :: weight
        character(160) :: error

        weight = sum(real(iatjpo(:nvec), real64))
        do kappa = minval(orbital_kappa(:orbital_count)),maxval(orbital_kappa(:orbital_count))
            ndim = count(orbital_kappa(:orbital_count) == kappa)
            if (ndim == 0) cycle
            allocate(orbits(ndim), matrix(ndim,ndim), rotation(ndim,ndim), occupation(ndim))
            orbits = pack([(i,i=1,orbital_count)], orbital_kappa(:orbital_count) == kappa)
            ! Enumerate actual orbitals, so gaps in n and n>20 are supported.
            do a = 2,ndim
                b = a
                do while (b > 1)
                    if (orbital_principal_n(orbits(b-1)) <= orbital_principal_n(orbits(b))) exit
                    temp = orbits(b-1)
                    orbits(b-1) = orbits(b)
                    orbits(b) = temp
                    b = b-1
                end do
            end do
            matrix = 0.0_real64
            do state = 1,nvec
                matrix = matrix+real(iatjpo(state),real64)*rho(orbits,orbits,state)/weight
            end do
            write(*,"(a,i5)") "KAPPA: ", kappa
            write(*,"(a,i5)") "NDIM: ", ndim
            write(*,*) "density matrix:"
            do a = 1,ndim
                write(*,"(i5,*(1x,es24.16))") orbital_principal_n(orbits(a)), matrix(:,a)
            end do
            call diagonalize_density(matrix, ordering, occupation, rotation, error)
            if (len_trim(error) /= 0) call density_fail(trim(error))
            points = maxval(orbital_last_grid_index(orbits))
            allocate(p(points,ndim), q(points,ndim), origin(ndim))
            p = radial_large_components(:points,orbits)
            q = radial_small_components(:points,orbits)
            origin = orbital_origin_parameters(orbits)
            call transform_radial_block(rotation, p, q, origin)
            radial_large_components(:points,orbits) = p
            radial_small_components(:points,orbits) = q
            orbital_origin_parameters(orbits) = origin
            orbital_last_grid_index(orbits) = points
            write(*,*) "occupations in output orbital order:"
            do a = 1,ndim
                write(*,"(i5,1x,es24.16)") orbital_principal_n(orbits(a)), occupation(a)
            end do
            write(*,*) "eigenvectors in output orbital order (columns printed as rows):"
            do a = 1,ndim
                write(*,"(i5,*(1x,es24.16))") orbital_principal_n(orbits(a)), rotation(:,a)
            end do
            deallocate(orbits, matrix, rotation, occupation, p, q, origin)
        end do
        open(newunit=unit, file=trim(name)//".nw", form="unformatted", status="replace", iostat=ios)
        if (ios /= 0) call density_fail("Cannot create natural-orbital file")
        write(unit) "G92RWF"
        do i = 1,orbital_count
            ! E is retained for GRASP file compatibility; it is not an NO energy.
            write(unit) orbital_principal_n(i), orbital_kappa(i), orbital_energy_parameters(i), orbital_last_grid_index(i)
            write(unit) orbital_origin_parameters(i), radial_large_components(:orbital_last_grid_index(i),i), &
               & radial_small_components(:orbital_last_grid_index(i),i)
            write(unit) radial_grid(:orbital_last_grid_index(i))
        end do
        close(unit, iostat=ios)
        if (ios /= 0) call density_fail("Error writing natural-orbital file")
        write(*,*) "Natural orbitals written to ", trim(name)//".nw"
    end subroutine make_natural_orbitals
end module density_calculation
