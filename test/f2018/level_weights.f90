program level_weights
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64
    use getaldwt_I, only: getaldwt
    implicit none (type, external)
    real(real64) :: guarded(5), expected(3)
    character(16) :: scenario
    call get_command_argument(1,scenario)
    select case(trim(scenario))
    case("equal")
        expected=1.0_real64/3
    case("standard")
        expected=[1.0_real64,3.0_real64,5.0_real64]/9
    case("custom","retry")
        expected=[0.2_real64,0.3_real64,0.5_real64]
    case default
        error stop "Unknown weight fixture"
    end select
    guarded=-99
    call getaldwt(3,guarded(2:4))
    if(.not.all(ieee_is_finite(guarded))) error stop "Nonfinite weights"
    if(maxval(abs(guarded(2:4)-expected))>1e-14_real64) error stop "Level weight values"
    if(abs(sum(guarded(2:4))-1.0_real64)>1e-14_real64) error stop "Level weight normalization"
    if(abs(guarded(1)+99.0_real64)>0.or.abs(guarded(5)+99.0_real64)>0) error stop "Weight extent"
    print *, "PASS level weights ",trim(scenario)
end program level_weights

! Controlled J labels: 2J+1 = 1,3,5. GETALDWT itself is the production body.
