program test_ranking
    use vast_kind_param, only: DOUBLE
    use iso_fortran_env, only: int64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use jj2lsj_postprocess, only: ls_ranking, rank_ls
    implicit none
    type(ls_ranking) :: ranked
    real(DOUBLE), allocatable :: values(:), old_weights(:), old_squares(:)
    integer, allocatable :: old_order(:)
    real(DOUBLE) :: old_norm, wa
    integer :: n, case_no, i, j, k
    character(len=16) :: argument
    call get_command_argument(1,argument)
    if (trim(argument)=='nonfinite') then
        call rank_ls([ieee_value(0.0_DOUBLE,ieee_quiet_nan)],ranked)
        stop 0
    end if
    do case_no=1,6
        n=257
        if (case_no==1) n=0
        if (case_no==2) n=1
        allocate(values(n),old_weights(n+1),old_squares(n+1),old_order(n+1))
        do i=1,n
            select case(case_no)
            case(2)
                values(i)=-0.0_DOUBLE
            case(3)
                values(i)=(-1.0_DOUBLE)**i*real(mod(i,7),DOUBLE)/8.0_DOUBLE
            case(4)
                values(i)=real(i,DOUBLE)/512.0_DOUBLE
            case(5)
                values(i)=real(n-i,DOUBLE)/512.0_DOUBLE
                if (mod(i,3)==0) values(i)=-0.0_DOUBLE
            case(6)
                values(i)=(-1.0_DOUBLE)**i/real(mod(37*i,263)+1,DOUBLE)
            end select
        end do
        ! Independent legacy insertion oracle, including its strict > tie rule.
        old_norm=0.0_DOUBLE
        old_weights=0.0_DOUBLE; old_squares=0.0_DOUBLE; old_order=0
        do i=1,n
            wa=values(i)
            old_norm=old_norm+wa*wa
            j=1
            do while (j<i)
                if (wa*wa>old_squares(j)) exit
                j=j+1
            end do
            do k=i,j+1,-1
                old_weights(k)=old_weights(k-1)
                old_squares(k)=old_squares(k-1)
                old_order(k)=old_order(k-1)
            end do
            old_weights(j)=wa; old_squares(j)=wa*wa; old_order(j)=i
        end do
        call rank_ls(values,ranked)
        if (any(ranked%order/=old_order(:n))) error stop 'unstable or incorrect ordering'
        if (transfer(ranked%norm,0_int64)/=transfer(old_norm,0_int64)) error stop 'norm order changed'
        do i=1,n
            if (transfer(values(ranked%order(i)),0_int64)/=transfer(old_weights(i),0_int64)) &
                error stop 'coefficient bits changed'
        end do
        deallocate(values,old_weights,old_squares,old_order)
    end do
    ! Large zero plateau: tied coefficients must retain their original LS IDs.
    allocate(values(84000))
    values=0.0_DOUBLE
    values(2)=-0.5_DOUBLE; values(84000)=0.5_DOUBLE
    call rank_ls(values,ranked)
    if (any(ranked%order(:3)/=[2,84000,1])) error stop 'large tie ordering'
    do i=4,84000
        if (ranked%order(i)/=i-1) error stop 'zero plateau order'
    end do
    print *, 'Stable ranking and original-order norms passed'
end program
