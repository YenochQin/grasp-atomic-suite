program level_weights
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use vast_kind_param, only: double
    use getaldwt_I, only: getaldwt
    implicit none
    real(double) :: guarded(5), expected(3)
    character(16) :: scenario
    call get_command_argument(1,scenario)
    select case(trim(scenario))
    case('equal')
        expected=1.0_double/3
    case('standard')
        expected=[1.0_double,3.0_double,5.0_double]/9
    case('custom','retry')
        expected=[0.2_double,0.3_double,0.5_double]
    case default
        error stop 'Unknown weight fixture'
    end select
    guarded=-99
    call getaldwt(3,guarded(2:4))
    if(.not.all(ieee_is_finite(guarded))) error stop 'Nonfinite weights'
    if(maxval(abs(guarded(2:4)-expected))>1e-14_double) error stop 'Level weight values'
    if(abs(sum(guarded(2:4))-1.0_double)>1e-14_double) error stop 'Level weight normalization'
    if(abs(guarded(1)+99.0_double)>0.or.abs(guarded(5)+99.0_double)>0) error stop 'Weight extent'
    print *, 'PASS level weights ',trim(scenario)
end program

! Controlled J labels: 2J+1 = 1,3,5. GETALDWT itself is the production body.
integer function itjpo(index)
    implicit none
    integer,intent(in) :: index
    if(index<1.or.index>3) error stop 'Unexpected level index'
    itjpo=2*index-1
end function
