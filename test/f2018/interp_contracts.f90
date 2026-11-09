program interp_contracts
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use vast_kind_param, only: double,long
    use interp_I, only: interp
    implicit none
    real(double) :: x(12),y(12),value,point,expected
    integer :: i,n
    character(16) :: scenario
    call get_command_argument(1,scenario)
    if(scenario=='legacy') then
        call legacy_parity()
        stop
    end if
    n=12
    do i=1,n
        x(i)=real(i,double)
    end do
    y=3+2*x+x*x
    select case(trim(scenario))
    case('regular')
        point=4.25_double;expected=29.5625_double
    case('left')
        point=1;expected=6
    case('short')
        n=2;y=2*x+1;point=1.5_double;expected=4
    case('zero')
        y=x-1.5_double;point=1.5_double;expected=0
    case('single','single_node')
        n=1;y=3;point=2;expected=3
        if(scenario=='single_node') point=1
    case('empty')
        n=0;point=0;expected=0
    case default
        error stop 'Unknown interpolation fixture'
    end select
    call interp(x(:n),y(:n),n,point,value,1e-12_double)
    if(scenario=='empty') error stop 'Empty interpolation table accepted'
    if(.not.ieee_is_finite(value)) error stop 'Nonfinite interpolation'
    if(abs(value-expected)>1e-12_double) error stop 'Interpolation analytic oracle'
    print '(A,Z16.16)', 'VALUE_BITS ',transfer(value,0_long)
    print *, 'PASS interpolation ',trim(scenario)
contains
    subroutine legacy_parity()
        real(double) :: nodes(40),values(40),points(0:40),point,result,expected
        integer :: sizes(3),n,i,j
        sizes=[11,17,40]
        do j=1,3
            n=sizes(j)
            do i=1,n
                nodes(i)=0.1_double*i+0.001_double*i*i
                values(i)=2+nodes(i)+nodes(i)**3/7
            end do
            points(0)=nodes(1)-0.02_double
            points(n)=nodes(n)+0.02_double
            do i=1,n-1
                points(i)=0.375_double*nodes(i)+0.625_double*nodes(i+1)
            end do
            do i=0,n
                point=points(i)
                call interp(nodes(:n),values(:n),n,point,result,1e-12_double)
                expected=2+point+point**3/7
                if(.not.ieee_is_finite(result)) error stop 'Nonfinite regular interpolation'
                if(abs(result-expected)>1e-12_double) error stop 'Regular interpolation polynomial'
                print '(2I4,1X,Z16.16)',n,i,transfer(result,0_long)
            end do
        end do
        print *, 'PASS interpolation legacy'
    end subroutine
end program
