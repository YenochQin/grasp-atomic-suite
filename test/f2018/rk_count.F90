


program rk_count
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: real64, int64
    use rk_fixture, only: calls
    use distributed_slater_integrals, only: local_slater_integral_indices, local_slater_integral_values, rk_local_count
    use orbital_csf_basis, only: orbital_count, orbital_l, orbital_two_j
    use tabulate_slater_integrals_interface, only: tabulate_slater_integrals
    implicit none (type, external)
    integer :: ranks, rank, count_total, jmax, total_calls, orbitals, i
    integer(int64) :: large_total
    nullify(local_slater_integral_indices,local_slater_integral_values)
    do orbitals=1,2
        orbital_count=orbitals
        orbital_l(1)=0; orbital_l(2)=1
        orbital_two_j(:orbital_count)=1
        do ranks=1,4
            total_calls=0
            do rank=0,ranks-1
                calls=0
                call tabulate_slater_integrals(rank,ranks,count_total,jmax)
                if(.not.all(ieee_is_finite(local_slater_integral_values))) error stop "Nonfinite integral cache"
                if(calls/=rk_local_count(count_total,rank,ranks)) error stop "Local Slater evaluation count"
                if(count(abs(local_slater_integral_values)>0.0_real64)/=calls) error stop "Local integral cache ownership"
                do i=1,count_total
                    if(mod(i,ranks)==rank) then
                        if(abs(local_slater_integral_values(i)-1.0_real64)>0) error stop "Owned integral"
                    else
                        if(abs(local_slater_integral_values(i))>0) error stop "Unowned integral"
                    end if
                end do
                total_calls=total_calls+calls
                deallocate(local_slater_integral_indices,local_slater_integral_values)
            end do
            if(total_calls/=count_total) error stop "Replicated cache counted more than once"
        end do
    end do
    large_total=0
    do rank=0,3
        large_total=large_total+int(rk_local_count(huge(0),rank,4),int64)
        if(rk_local_count(0,rank,4)/=0) error stop "Empty integral cache"
    end do
    if(large_total/=int(huge(0),int64)) error stop "Count arithmetic overflow"
    print *, "PASS production integral partition, empty ranks, repeated calls and count boundary"
end program rk_count
