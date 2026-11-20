module evaluate_breit_hartree_z_interface
    implicit none (type, external)

    ! Module interface
    private
    public :: evaluate_breit_hartree_z
contains

    !***********************************************************************
    !                                                                      *
    subroutine evaluate_breit_hartree_z(k, i, j)
        ! Module dependencies
        !                                                                      *
        !   This subroutine evaluates Hartree Z-functionals:                   *
        !                                                                      *
        !              (k)                     k                               *
        !             Z   [f(r);r] =  I ( (s/r)   f(s) ; 0 - r )               *
        !                                                                      *
        !   where  I ( g(r,s) ; range )  denotes the integral of g(r,s) over   *
        !   range  in  s .    The Z-functional is tabulated in  COMMON/TATB/   *
        !   in array  TB . The f-function is assumed tabulted in array  TA .   *
        !                                                                      *
        !   Written by Farid A Parpia, at Oxford   Last updated: 14 Oct 1992   *
        !                                                                      *
        !***********************************************************************
        use, intrinsic :: iso_fortran_env, only: real64
        use grasp_parameters,   only: radial_workspace_capacity
        use radial_difference_coefficients, only: five_point_weights
        use radial_grid_state, only: radial_point_count, radial_grid, radial_grid_derivative
        use radial_quadrature_weights, only: endpoint_weight, odd_interior_weight, middle_weight
        use radial_quadrature_workspace, only: radial_workspace_b, integration_last_index, &
            radial_workspace_capacity, radial_workspace_a
        implicit none (type, external)

        ! Arguments
        integer, intent(in) :: k
        integer  :: i
        integer  :: j    !!!! Arument not referenced !!!

        ! Local declarations
        integer :: ii, mtpp1, mtpp3, mtpp4, kk
        real(real64) :: rhop(radial_workspace_capacity), rttk(radial_workspace_capacity), temp(radial_workspace_capacity)
        real(real64) :: sum, zklim
        !-----------------------------------------------
        !
        if (k > 0) then
            rttk(2:radial_point_count) = radial_grid(2:radial_point_count)**k
        end if
        !
        !   MTP is fed in through  COMMON/TATB/
        !
        mtpp1 = integration_last_index + 1
        mtpp3 = integration_last_index + 3
        mtpp4 = integration_last_index + 4
        !
        !   Compute  RP(S)*F(S)  and store it in  RHOP
        !
        rhop(2:integration_last_index) = radial_grid_derivative(2:integration_last_index)* &
            radial_workspace_a(2:integration_last_index)
        !
        !   Fill array TEMP with r**k * RHOP
        !
        temp(1) = 0.0_real64
        if (k == 0) then
            temp(2:integration_last_index) = rhop(2:integration_last_index)
        else
            temp(2:integration_last_index) = rttk(2:integration_last_index)*rhop(2:integration_last_index)
        end if
        !
        !   Set an additional four points to zero
        !
        temp(mtpp1:mtpp4) = 0.0_real64
        !
        !                                     k
        !   Compute the first few values of  r  * ZK  using semi-open
        !   Newton-Cotes formulae. Store the Breit Z^k in workspace B.
        !
        radial_workspace_b(1) = 0.0_real64
        do ii = 2, 4
            sum = 0.0_real64
            do kk = 2, 5
                sum = sum + five_point_weights(kk,ii)*temp(kk)
            end do
            radial_workspace_b(ii) = sum
        end do
        !                         k
        !   Compute remainder of r  * ZK: march out to MTP+3
        !
        do ii = 5, mtpp3
            radial_workspace_b(ii) = radial_workspace_b(ii-4) + endpoint_weight*(temp(ii-4)+temp(ii)) + &
                odd_interior_weight*(temp(ii-3)+temp(ii-1)) + middle_weight*temp(ii-2)
        end do
        !                                       k   (k)
        !   Determine the asymptotic value of  r * Z
        !
        !   Compute ZK
        !
        zklim = radial_workspace_b(mtpp3)
        !
        if (k == 0) then
            !
            radial_workspace_b(mtpp4:radial_point_count) = zklim
        !
        else
            !
            radial_workspace_b(2:mtpp3) = radial_workspace_b(2:mtpp3)/rttk(2:mtpp3)
            !
            radial_workspace_b(mtpp4:radial_point_count) = zklim/rttk(mtpp4:radial_point_count)
            !
        end if
        !
        return
    end subroutine evaluate_breit_hartree_z
end module evaluate_breit_hartree_z_interface
