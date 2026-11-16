program bndinv_precision
    use iso_fortran_env, only: real32, real64
    use bndinv_i
    implicit none (type, external)
    real(real64) :: a(2,2), work(2,2), original(2,2), determinant, values(5), x
    integer :: i, status
    values = [nearest(1e-8_real64,-1.0_real64), 1e-8_real64, &
              nearest(1e-8_real64,1.0_real64), &
              (real(1e-8_real32,real64)+1e-8_real64)/2, 0.25_real64]
    do i=1,size(values)
        x=values(i)
        a=0
        a(1,1)=1
        a(2,2)=x
        call bndinv(a,work,2,determinant,0.0_real64,status,2)
        if (status /= merge(1,0,x < 1e-8_real64)) error stop "BNDINV default boundary"
        if (status == 0) then
            if (abs(a(2,2)*x-1.0_real64) > 4*epsilon(x)) error stop "BNDINV diagonal inverse"
        end if
        a=0
        a(1,1)=1
        a(2,2)=x
        call bndinv(a,work,2,determinant,0.5_real64,status,2)
        if (status /= 1) error stop "BNDINV explicit threshold ignored"
    end do
    original=reshape([4.0_real64,1.0_real64,1.0_real64,3.0_real64],[2,2])
    a=original
    call bndinv(a,work,2,determinant,0.0_real64,status,2)
    if (status /= 0) error stop "BNDINV normal matrix rejected"
    if (maxval(abs(matmul(a,original)-reshape([1,0,0,1],[2,2]))) > &
        32*epsilon(1.0_real64)) error stop "BNDINV inverse residual"
    print *, "BNDINV: boundary, explicit tolerance and normal inverse passed"
end program bndinv_precision
