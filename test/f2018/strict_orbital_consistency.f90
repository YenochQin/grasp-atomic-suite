! Real strict IMPROV/CONSIS/DAMPOR and QUAD; only radial equation assembly
! and SOLVE are controlled. This catches the comparison's placement in IMPROV.
module orbital_fixture
    use vast_kind_param, only: double
    use parameter_def, only: nnnp
    implicit none
    real(double) :: reference(nnnp), proposal_scale, angle
    integer :: damping_calls=0
end module

program strict_orbital_consistency
    use orbital_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use grid_C, only: n, h, r, rp
    use def_C, only: accy, nsic
    use wave_C, only: pf, qf, pz, mf
    use scf_C, only: ucf, method, scnsty
    use orb_C, only: e, gama, nw, np, nh
    use damp_C, only: odamp
    use invt_C, only: noinvt
    use orthct_C, only: orthst
    use mpi_C, only: myid, nprocs
    use improvmpi_I, only: improvmpi
    use setqic_I, only: setqic
    use rint_I, only: rint
    implicit none
    integer :: i, scenario, orbital
    real(double) :: damping_max, expected, saved_change
    logical :: eol, sorting

    ! P=r/sqrt(243), Q=0 on [0,9]: analytic norm integral r^2/243 = 1.
    ! The production QUAD is exact for this quadratic on complete panels.
    n=10
    h=1.0_double
    r=0
    rp=1
    do i=1,n
        r(i)=real(i-1,double)
    end do
    call setqic()
    reference=0
    reference(:n)=r(:n)/sqrt(243.0_double)
    allocate(pf(nnnp,1),qf(nnnp,1))
    nw=1
    orbital=1
    np(1)=1
    nh(1)='1s'
    gama(1)=1
    e(1)=1
    ucf(1)=2
    method(1)=3
    accy=1.5625e-8_double
    nsic=0
    orthst=.false.
    noinvt(1)=.true.
    myid=0
    nprocs=1
    eol=.true.
    sorting=.false.
    do scenario=1,4
        pf(:,1)=reference
        qf=0
        pz(1)=1.0_double/sqrt(243.0_double)
        mf(1)=n
        if (abs(rint(orbital,orbital,0)-1.0_double)>1e-14_double) error stop 'Fixture norm'
        proposal_scale=2.0_double
        angle=0
        if (scenario==2) proposal_scale=0.25_double
        if (scenario==3) angle=1e-4_double
        if (scenario==4) angle=accy/4
        ! A small accepted step must not conceal a larger undamped proposal.
        odamp(1)=0.9999_double
        damping_calls=0
        damping_max=0
        call improvmpi(eol,orbital,sorting,damping_max)
        expected=sqrt(ucf(1))*maxval(reference)*(abs(cos(angle)-1.0_double)+abs(sin(angle)))
        if (.not.ieee_is_finite(scnsty(1))) error stop 'Nonfinite consistency'
        if (abs(scnsty(1)-expected)>1e-13_double) error stop 'Normalized proposal consistency'
        if (abs(rint(orbital,orbital,0)-1.0_double)>1e-14_double) error stop 'Updated orbital norm'
        saved_change=sqrt(ucf(1))*maxval(abs(pf(:,1)-reference)+abs(qf(:,1)))
        if (scenario==3) then
            if (scnsty(1)<=accy.or.damping_calls/=1) error stop 'Shape change incorrectly accepted'
            if (saved_change>=accy) error stop 'Fixture must hide shape change after damping'
        else
            if (scnsty(1)>accy.or.damping_calls/=0) error stop 'Scale/below-threshold change rejected'
        end if
    end do
    print *, 'PASS normalized undamped proposal; scale invariance and shape rejection'
end program

subroutine solve(j,fail,inv,jp,nnp)
    use orbital_fixture
    use int_C, only: p,q,p0,mtp0
    use grid_C, only: n
    implicit none
    integer :: j
    logical :: fail
    integer, intent(inout) :: inv
    integer, intent(out) :: jp,nnp
    if(j/=1.or.inv/=0) error stop 'SOLVE fixture arguments'
    p=proposal_scale*cos(angle)*reference
    q=proposal_scale*sin(angle)*reference
    p0=proposal_scale*cos(angle)/sqrt(243.0_double)
    mtp0=n
    jp=4
    nnp=0
    fail=.false.
end subroutine

subroutine cofpotmpi(eol,j,npts)
    implicit none
    logical :: eol
    integer :: j,npts
    if(.not.eol.or.j/=1.or.npts/=10) error stop 'COFPOT fixture arguments'
end subroutine

subroutine defcor(j)
    implicit none
    integer,intent(in) :: j
    if(j/=1) error stop 'DEFCOR fixture argument'
end subroutine

subroutine dampck(ipr,j,ed1,ed2)
    use orbital_fixture, only: double,damping_calls
    implicit none
    integer,intent(inout) :: ipr
    integer,intent(in) :: j
    real(double),intent(inout) :: ed1,ed2
    if(j/=1.or.ipr/=0.or.abs(ed1)>0.or.abs(ed2-1.0_double)>0) error stop 'DAMPCK fixture arguments'
    damping_calls=damping_calls+1
end subroutine

subroutine setlagmpi(eol)
    implicit none
    logical :: eol
    if(eol) error stop 'Unexpected SOLVE retry'
    error stop 'Unexpected SETLAG'
end subroutine

subroutine orthy(nw,jp,lsort)
    implicit none
    integer,intent(in) :: nw,jp
    logical,intent(in) :: lsort
    if(nw/=1.or.jp/=1.or.lsort) error stop 'ORTHY fixture arguments'
    error stop 'Unexpected ORTHY'
end subroutine
