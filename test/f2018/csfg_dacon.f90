program csfg_dacon
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use vast_kind_param, only: double
    use parameter_def, only: nnnp
    use debug_C, only: ldbpr
    use def_C, only: c
    use grid_C, only: h,r,rpor
    use npot_C, only: zz
    use orb_C, only: np,nh,nak
    use pote_C, only: xp,xq
    use csfg_scf_C, only: ndcof,nda,da,ldbpdacon
    use wave_C, only: pf,qf,mf
    use dacon_I, only: dacon
    implicit none
    integer :: i, mode
    real(double) :: reference_p(5),reference_q(5)
    if(any(ldbpdacon)) error stop 'DACON diagnostics must default off'
    ldbpr=.false.
    c=2; h=0.5_double; r=1; rpor=2; zz=2
    np(1)=1; nh(1)='1s'; nak(1)=-1; mf(1)=5
    allocate(pf(nnnp,1),qf(nnnp,1),nda(1),da(1))
    pf=0; qf=0.5_double
    do i=1,5
        pf(i,1)=real(i,double)
    end do
    nda=1;da=0.5_double;ndcof=1
    open(99,status='scratch',form='formatted')
    do mode=0,3
        ldbpdacon=[btest(mode,0),btest(mode,1)]
        xp=7;xq=-3
        call dacon(1)
        if(.not.all(ieee_is_finite(xp(:5))).or..not.all(ieee_is_finite(xq(:5)))) error stop 'Nonfinite DACON'
        do i=2,5
            if(abs(xp(i)-(7.25_double-0.5_double*i))>1e-14_double) error stop 'DACON P oracle'
            if(abs(xq(i)-(-3.75_double+0.5_double*i))>1e-14_double) error stop 'DACON Q oracle'
        end do
        if(abs(xp(1)-7.0_double)>0.or.abs(xq(1)+3.0_double)>0) error stop 'DACON origin'
        if(mode==0) then
            reference_p=xp(:5);reference_q=xq(:5)
        else
            if(maxval(abs(xp(:5)-reference_p))+maxval(abs(xq(:5)-reference_q))>0) &
                error stop 'DACON diagnostics changed numerical values'
        end if
    end do
    close(99)
    print *, 'PASS production CSFG DACON'
end program

subroutine dpbdt(j)
    use tatb_C, only: ta,tb
    implicit none
    integer,intent(in) :: j
    if(j/=1) error stop 'DACON derivative fixture'
    ta=3;tb=-2
end subroutine
