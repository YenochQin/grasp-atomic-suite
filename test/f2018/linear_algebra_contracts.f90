program linear_algebra_contracts
    use, intrinsic :: iso_fortran_env, only: real32
    use, intrinsic :: iso_fortran_env, only: real64
    use blas_lapack_interfaces, only: dcopy, ddot, dspevx, slamch, dlamch
    implicit none (type, external)
    real(real64), target :: backing(8)
    real(real64), pointer :: view(:)
    real(real64) :: copied(3), x(3), y(3), ap(3), w(2), z(2,2), work(16), h(2,2)
    integer :: found, info, iwork(10), ifail(2), i

    if (storage_size(0) /= 32 .or. storage_size(0.0_real64) /= 64) error stop "LP64 contract"
    if (slamch("S") /= tiny(0.0_real32)) error stop "SLAMCH return ABI"
    if (dlamch("S") /= tiny(0.0_real64)) error stop "DLAMCH return ABI"
    x = [1.25_real64,2.0_real64,3.0_real64]
    y = [-2.0_real64,4.0_real64,1.25_real64]
    if (ddot(3,x,1,y,1) /= 9.25_real64) error stop "DDOT return ABI"
    backing = [0.0_real64,1.0_real64,2.0_real64,3.0_real64,4.0_real64,5.0_real64,6.0_real64,7.0_real64]
    view => backing
    call dcopy(3,view(3:),1,copied,1)
    if (any(copied /= [2.0_real64,3.0_real64,4.0_real64])) error stop "pointer subarray copy"
    call dcopy(3,view(3:),-1,copied,1)
    if (any(copied /= [4.0_real64,3.0_real64,2.0_real64])) error stop "negative stride copy"
    view => backing(::2)
    call dcopy(3,view,1,copied,1)
    if (any(copied /= [0.0_real64,2.0_real64,4.0_real64])) error stop "strided view copy"

    h = reshape([2.0_real64,1.0_real64,1.0_real64,2.0_real64],[2,2])
    ap = [2.0_real64,1.0_real64,2.0_real64]
    call dspevx("V","A","U",2,ap,0.0_real64,0.0_real64,1,2,0.0_real64,found,w,z,2,work,iwork,ifail,info)
    if (info /= 0 .or. found /= 2) error stop "DSPEVX failed"
    if (maxval(abs(w-[1.0_real64,3.0_real64])) > 1e-12_real64) error stop "DSPEVX eigenvalues"
    do i=1,2
        if (maxval(abs(matmul(h,z(:,i))-w(i)*z(:,i))) > 1e-12_real64) error stop "DSPEVX residual"
    end do
    if (abs(dot_product(z(:,1),z(:,2))) > 1e-12_real64) error stop "DSPEVX orthogonality"
end program linear_algebra_contracts
