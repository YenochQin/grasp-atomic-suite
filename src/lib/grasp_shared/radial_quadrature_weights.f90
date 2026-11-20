module radial_quadrature_weights
    ! Module dependencies
    use, intrinsic :: iso_fortran_env, only: real64
    implicit none (type, external)

    ! Module interface
    private :: real64

    ! Module state
    ! One five-point panel uses (c1,c2,c3,c2,c1). Adjacent panels share an endpoint of weight c4.
    real(real64) :: endpoint_weight, odd_interior_weight, middle_weight, shared_endpoint_weight
end module radial_quadrature_weights
