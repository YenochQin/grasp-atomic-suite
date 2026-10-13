program start_trace_parity
    use vast_kind_param, only: double
    use parameter_def, only: nnnp
    use orb_C, only: np, nh, nak, e, gama
    use grid_C, only: r, rp, h
    use def_C, only: c, z, accy
    use npar_C, only: nparm
    use pote_C, only: xp, xq, yp
    use wave_C, only: pf, qf, pz
    use scf_C, only: ndcof, nda, da
    use cnc_C, only: cnc6c
    use p0_trace, only: p0_trace_cycle, p0_trace_update, p0_trace_iteration
    use start_trace
    use start_I
    use start_diag_I
    implicit none
    real(double) :: original_p(nnnp), original_q(nnnp), traced_p(nnnp), traced_q(nnnp)
    real(double) :: original_q0, traced_q0
    integer :: rank, mode, nucleus, orbital, itype, i
    character(len=1024) :: directory
    character(len=16) :: argument
    call get_command_argument(1, directory)
    call get_command_argument(2, argument)
    read(argument,*) rank
    allocate(pf(nnnp,2), qf(nnnp,2), nda(1), da(1))
    pf = 0.01d0
    qf = -0.001d0
    pz(:2) = [221.d0, 0.2d0]
    np(:2) = 3
    nh(:2) = ['d ', 'd-']
    nak(:2) = [-3,2]
    gama(:2) = [3.d0,2.d0]
    e(:2) = [2.d0,3.d0]
    ndcof = 1
    nda(1) = 1
    da(1) = 0.02d0
    c = 137.035999139d0
    z = 28.d0
    h = 0.05d0
    accy = h**6
    do i = 1, nnnp
        r(i) = 1.d-6*(exp((i-1)*h)-1.d0)
        rp(i) = r(i) + 1.d-6
    end do
    xp = 1.d-22
    xq = -2.d-23
    yp = 0.01d0
    call start_trace_open(trim(directory), rank)
    do mode = 1, 2
        cnc6c = 0.d0
        if (mode == 2) cnc6c = 100.d0 ! Force the inherited 36-iteration fallback.
        call p0_trace_cycle(mode)
        do nucleus = 0, 2, 2
            nparm = nucleus
            do orbital = 1, 2
                call p0_trace_update(orbital)
                do itype = 1, 3
                    call p0_trace_iteration(itype)
                    original_p = 0.d0
                    original_q = 0.d0
                    traced_p = original_p
                    traced_q = original_q
                    call start(orbital, itype, pz(orbital), original_p, original_q0, original_q)
                    call start_diag(orbital, itype, pz(orbital), traced_p, traced_q0, traced_q)
                    if (any(original_p /= traced_p) .or. any(original_q /= traced_q)) &
                        error stop 'START diagnostic changed a radial value'
                    if (original_q0 /= traced_q0) error stop 'START diagnostic changed Q0'
                end do
            end do
        end do
    end do
    call start_trace_close()
    print *, 'START parity checks passed.'
end program start_trace_parity
