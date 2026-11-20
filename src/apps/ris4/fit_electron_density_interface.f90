module fit_electron_density_interface
    implicit none (type, external)

    ! Procedure interfaces
    interface
        !...Translated by Gediminas Gaigalas 11/18/19
        subroutine fit_electron_density (xvec,yvec,z,par,nrnuc,f,rho,res)
            ! Module dependencies
            use, intrinsic :: iso_fortran_env, only: real64
            use grasp_parameters,    only: max_radial_points, radial_workspace_capacity
            implicit none (type, external)

            ! Arguments and local declarations
            real(real64) :: xvec(radial_workspace_capacity)
            real(real64) :: yvec(max_radial_points)
            !GG      REAL(DOUBLE), DIMENSION(NNNP) :: XVEC, YVEC
            real(real64) :: p(5), f(5)
            real(real64) :: par(2)
            real(real64) :: z, drms, rho, res
            integer      :: nrnuc, nparfit
        end subroutine fit_electron_density
    end interface
end module fit_electron_density_interface
