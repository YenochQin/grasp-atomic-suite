!***********************************************************************
!                                                                      *
real(real64) function fullerton_rinker_k_function (x, n)
    ! Module dependencies
    !                                                                      *
    !   This  function  evaluates the KN(X) functions using the analytic   *
    !   functions defined in tables 1 and 3 of Fullerton and Rinker.       *
    !                                                                      *
    !                                           Last update: 09 Oct 1992   *
    !                                                                      *
    !***********************************************************************
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Arguments
    integer, intent(in) :: n
    real(real64), intent(in) :: x

    ! Local declarations
    integer :: np(4)
    integer :: k, nn, i
    real(real64) :: p(10, 4)
    real(real64) :: b(2, 4)
    real(real64) :: c(3, 4)
    real(real64) :: d(5, 4), e(5, 4)
    real(real64) :: xn, sum, x2, bsum, csum, dsum, esum, xm

    ! Initial values
    !-----------------------------------------------
    !
    !
    data p/ 8.8357293375e-1_real64, -2.8259817381e-1_real64, -5.8904879578e-1_real64, &
        1.2500133434e-1_real64, -3.2729913852e-2_real64, 8.2888574511e-3_real64, &
        -1.0327765800e-5_real64, 6.3643668900e-5_real64, 0.0_real64, 0.0_real64, -7.1740181754e-1_real64&
        , 1.1780972274_real64, -3.7499963087e-1_real64, 1.3089675530e-1_real64, &
        -3.8258286439e-2_real64, -2.4297287300e-5_real64, -3.5920148670e-4_real64, &
        -1.7170090700e-5_real64, 0.0_real64, 0.0_real64, 9.9999999987e-1_real64, 1.9770200000e-8_real64&
        , -7.5000050190e-1_real64, 7.8540306316e-1_real64, -3.4988601655e-1_real64, &
        6.4596333000e-5_real64, -9.8189080747e-3_real64, 8.6513145800e-5_real64, &
        -2.3969236620e-4_real64, 0.0_real64, 6.0000000002_real64, -6.4305200000e-8_real64, &
        2.1049413000e-6_real64, -2.6711271500e-5_real64, -1.3705236152e-1_real64, &
        -6.3476104090e-4_real64, -7.8739801501e-2_real64, -1.9641740173e-3_real64, &
        -3.4752369349e-3_real64,  - 7.3145316220e-4_real64/
    !
    data b/ -3.19999594323e2_real64, 2.53900995981_real64, -6.40514843293e1_real64, &
        7.11722714285e-1_real64, 5.19010136460e3_real64, 8.28495496200e1_real64, &
        3.18150793824e2_real64, 4.33898867347e1_real64/
    !
    data c/ -3.19999594333e2_real64, 2.53901020662_real64, 0.0_real64, 6.40514843287e1_real64&
        , -7.11722686403e-1_real64, 8.042207748e-4_real64, 2.76805406060e4_real64, &
        -3.27039477790e2_real64, 0.0_real64, 8.48402116837e2_real64, -2.56939867765e1_real64, &
        3.20844906346e-1_real64/
    !
    data d/ 5.018065179_real64, 7.151891262e1_real64, 2.116209929e2_real64, &
        3.140327478e1_real64, -1.0_real64, 2.172386409e2_real64, 1.643364528e3_real64, &
        2.122244512e3_real64, -4.512004044e1_real64, 1.0_real64, 8.540770444_real64, &
        6.076242766e1_real64, 9.714630584e1_real64, 3.154973593e1_real64, 1.0_real64, &
        5.9243015865e-1_real64, 2.0596312871_real64, 3.7785190424_real64, 3.5614853214_real64, &
        1.0_real64/
    !
    data e/ 2.669207401_real64, 5.172549669e1_real64, 2.969809720e2_real64, &
        5.364324164e2_real64, 1.535335924e2_real64, 1.155589983e2_real64, 1.292191441e3_real64, &
        3.831198012e3_real64, 2.904410075e3_real64, 0.0_real64, 4.543392478_real64, &
        3.514920169e1_real64, 6.019668656e1_real64, 8.468839579_real64, 0.0_real64, &
        3.1511867816e-1_real64, 3.473245222e-1_real64, 3.8791936870e-2_real64, -1.3059741497e-3_real64&
        , 0.0_real64/
    !
    data np/ 8, 8, 9, 10/
    !
    if (x /= 0.0_real64) then
        select case (n)
        case (0, 1)
            k = n + 1
            xn = 1.0_real64
        case (3)
            k = n
            xn = 1.0_real64/x**2
        case (5)
            k = n - 1
            xn = 1.0_real64/x**4
        case default
            write (*, '(/," Attempt to calculate FUNK (X,N) for N other than"," 0, 1, 3 and 5.")')
            stop
        end select
        if (x <= 1.0_real64) then
            !
            !   Calculate function for X < = 1
            !
            nn = np(k)
            sum = 0.0_real64
            do i = 1, nn
                sum = sum + p(i,k)*xn
                xn = xn*x
            end do
            x2 = x*x
            bsum = b(1,k) + x2*(b(2,k)+x2)
            csum = c(1,k) + x2*(c(2,k)+x2*c(3,k))
            if (k == 1) bsum = bsum*x
            if (k == 3) bsum = bsum*x2
            sum = sum + bsum*dlog(x)/csum
            fullerton_rinker_k_function = sum
            return
        end if
        !
        !   Calculate function for X > 1
        !
        xn = 1.0_real64
        dsum = 0.0_real64
        esum = 0.0_real64
        do i = 1, 5
            dsum = dsum + d(i,k)*xn
            esum = esum + e(i,k)*xn
            xn = xn/x
        end do
        xm = -x
        sum = dsum*exp(xm)/(esum*dsqrt(x**3))
        fullerton_rinker_k_function = sum
        return
    end if
    if (n /= 0) then
        write (*, '(/," Attempt to calculate FUNK (0,N) for N > 0")')
        stop
    end if
    fullerton_rinker_k_function = p(1,1)
end function fullerton_rinker_k_function
