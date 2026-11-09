module rk_fixture
    implicit none
    integer :: calls=0
end module

program rk_count
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use vast_kind_param, only: double,long
    use rk_fixture
    use cteilsrk_C, only: indteirk,valteirk,rk_local_count
    use orb_C, only: nw,nkl,nkj
    use genintrk_I, only: genintrk
#ifdef RK_CSFG
    use symmatrix_mod, only: kstartia
#endif
    implicit none
    integer :: ranks, rank, count_total, jmax, total_calls, orbitals, i
    integer(long) :: large_total
    nullify(indteirk,valteirk)
    do orbitals=1,2
        nw=orbitals
        nkl(1)=0; nkl(2)=1
        nkj(:nw)=1
        do ranks=1,4
            total_calls=0
            do rank=0,ranks-1
                calls=0
                call genintrk(rank,ranks,count_total,jmax)
                if(.not.all(ieee_is_finite(valteirk))) error stop 'Nonfinite integral cache'
                if(calls/=rk_local_count(count_total,rank,ranks)) error stop 'Local Slater evaluation count'
#ifdef RK_CSFG
                if(kstartia(jmax*nw+nw)/=count_total) error stop 'CSFG integral cache boundary'
#endif
                if(count(abs(valteirk)>0.0_double)/=calls) error stop 'Local integral cache ownership'
                do i=1,count_total
                    if(mod(i,ranks)==rank) then
                        if(abs(valteirk(i)-1.0_double)>0) error stop 'Owned integral'
                    else
                        if(abs(valteirk(i))>0) error stop 'Unowned integral'
                    end if
                end do
                total_calls=total_calls+calls
                deallocate(indteirk,valteirk)
            end do
            if(total_calls/=count_total) error stop 'Replicated cache counted more than once'
        end do
    end do
    large_total=0
    do rank=0,3
        large_total=large_total+int(rk_local_count(huge(0),rank,4),long)
        if(rk_local_count(0,rank,4)/=0) error stop 'Empty integral cache'
    end do
    if(large_total/=int(huge(0),long)) error stop 'Count arithmetic overflow'
    print *, 'PASS production integral partition, empty ranks, repeated calls and count boundary'
end program

function slater(ia,ib,ic,id,k) result(value)
    use vast_kind_param, only: double
    use rk_fixture, only: calls
    implicit none
    integer,intent(in) :: ia,ib,ic,id,k
    real(double) :: value
    if(min(ia,ib,ic,id)<1.or.k<0) error stop 'Invalid Slater fixture indices'
    calls=calls+1
    value=1.0_double
end function
