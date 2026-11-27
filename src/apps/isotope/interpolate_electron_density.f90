! Interpolate the electron density used for field-shift factors.
module interpolate_electron_density_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: interpolate_electron_density
contains

    !   This routine uses interpolating polynomial d(r) = d0 + d2*r^2 + d4*r^4
    !   to extrapolate electron density at r = 0.
    !
    ! Written by Jorgen Ekman Jul. 2016
    subroutine interpolate_electron_density(xa,ya,dens)
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        implicit none (type, external)

        ! Arguments
        real(real64), intent(out) :: dens
        real(real64),  intent(in) :: xa(3), ya(3)

        ! Local declarations
        real(real64) :: m(3, 3), mi(3, 3)
        real(real64) :: rm(3), pm(3)

        real(real64) :: mdet

        rm(1) = ya(1)
        rm(2) = ya(2)
        rm(3) = ya(3)

        m(1,1) = 1.0_real64
        m(1,2) = xa(1)**2.0_real64
        m(1,3) = xa(1)**4.0_real64
        m(2,1) = 1.0_real64
        m(2,2) = xa(2)**2.0_real64
        m(2,3) = xa(2)**4.0_real64
        m(3,1) = 1.0_real64
        m(3,2) = xa(3)**2.0_real64
        m(3,3) = xa(3)**4.0_real64
        mdet = m(1,1)*m(2,2)*m(3,3)-m(1,1)*m(2,3)*m(3,2)- &
            m(1,2)*m(2,1)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+ &
            m(1,3)*m(2,1)*m(3,2)-m(1,3)*m(2,2)*m(3,1)

        ! DETERMINE INVERSE OF MATRIX
        mi(1,1) = (m(2,2)*m(3,3)-m(2,3)*m(3,2))/mdet
        mi(1,2) = (m(1,3)*m(3,2)-m(1,2)*m(3,3))/mdet
        mi(1,3) = (m(1,2)*m(2,3)-m(1,3)*m(2,2))/mdet

        mi(2,1) = (m(2,3)*m(3,1)-m(2,1)*m(3,3))/mdet
        mi(2,2) = (m(1,1)*m(3,3)-m(1,3)*m(3,1))/mdet
        mi(2,3) = (m(1,3)*m(2,1)-m(1,1)*m(2,3))/mdet

        mi(3,1) = (m(2,1)*m(3,2)-m(2,2)*m(3,1))/mdet
        mi(3,2) = (m(1,2)*m(3,1)-m(1,1)*m(3,2))/mdet
        mi(3,3) = (m(1,1)*m(2,2)-m(1,2)*m(2,1))/mdet

        ! DETERMINE PARAMETERS FROM FIT
        pm(1) = mi(1,1)*rm(1)+mi(1,2)*rm(2)+mi(1,3)*rm(3)
        pm(2) = mi(2,1)*rm(1)+mi(2,2)*rm(2)+mi(2,3)*rm(3)
        pm(3) = mi(3,1)*rm(1)+mi(3,2)*rm(2)+mi(3,3)*rm(3)

        ! Finally RHO(0) in au^{-3} is defined
        dens = pm(1)

        return
    end subroutine interpolate_electron_density
end module interpolate_electron_density_interface
