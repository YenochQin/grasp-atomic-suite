! Initialize the derivative and Newton-Cotes weights used by radial integration and differentiation.
module radial_stencil_weights
    implicit none (type, external)

    ! Module interface
    private
    public :: setqic
contains
    subroutine setqic
        ! Module dependencies
        use, intrinsic :: iso_fortran_env, only: real64
        use radial_difference_coefficients, only: five_point_weights => cnc5c, six_point_weights => cnc6c
        use radial_grid_state, only: radial_workspace_capacity, radial_point_count, radial_grid, radial_grid_derivative, &
            radial_derivative_over_radius, radial_grid_scale, transformed_grid_step, asymptotic_grid_step
        use radial_derivative_stencil, only: derivative_weights => a
        use radial_quadrature_weights, only: endpoint_weight => c1, odd_interior_weight => c2, middle_weight => c3, &
            & shared_endpoint_weight => c4
        use radial_stepping_coefficients, only: starting_step_weights => c
        implicit none (type, external)

        ! Local declarations
        integer :: weight_index, column_index
        real(real64) :: derivative_numerators(13, 13)
        real(real64) :: starting_weights(6)
        real(real64) :: five_point_numerators(5, 2:5)
        real(real64) :: six_point_numerators(6, 2:6)
        real(real64) :: derivative_denominator, starting_denominator, five_point_denominator, six_point_denominator, weight_scale
        logical :: initialize_derivative_table

        ! Initial values
        data (derivative_numerators(1,weight_index),weight_index=1,13)/ -1486442880.0_real64, 5748019200.0_real64, &
            -15807052800.0_real64, 35126784000.0_real64, -59276448000.0_real64, &
            75873853440.0_real64, -73766246400.0_real64, 54195609600.0_real64, &
            -29638224000.0_real64, 11708928000.0_real64, -3161410560.0_real64, 522547200.0_real64&
            ,  - 39916800.0_real64/
        data (derivative_numerators(2,weight_index),weight_index=1,13)/ -39916800.0_real64, -967524480.0_real64, &
            & 2634508800.0_real64&
            , -4390848000.0_real64, 6586272000.0_real64, -7903526400.0_real64, 7376624640.0_real64&
            , -5269017600.0_real64, 2822688000.0_real64, -1097712000.0_real64, 292723200.0_real64&
            , -47900160.0_real64, 3628800.0_real64/
        data (derivative_numerators(3,weight_index),weight_index=1,13)/ 3628800.0_real64, -87091200.0_real64, &
            & -684478080.0_real64, &
            1596672000.0_real64, -1796256000.0_real64, 1916006400.0_real64, -1676505600.0_real64, &
            1149603840.0_real64, -598752000.0_real64, 228096000.0_real64, -59875200.0_real64, &
            9676800.0_real64,  - 725760.0_real64/
        data (derivative_numerators(4,weight_index),weight_index=1,13)/ -725760.0_real64, 13063680.0_real64, &
            & -143700480.0_real64, &
            -476910720.0_real64, 1077753600.0_real64, -862202880.0_real64, 670602240.0_real64, &
            -431101440.0_real64, 215550720.0_real64, -79833600.0_real64, 20528640.0_real64, &
            -3265920.0_real64, 241920.0_real64/
        data (derivative_numerators(5,weight_index),weight_index=1,13)/ 241920.0_real64, -3870720.0_real64, 31933440.0_real64, &
            -212889600.0_real64, -303937920.0_real64, 766402560.0_real64, -447068160.0_real64, &
            255467520.0_real64, -119750400.0_real64, 42577920.0_real64, -10644480.0_real64, &
            1658880.0_real64,  - 120960.0_real64/
        data (derivative_numerators(6,weight_index),weight_index=1,13)/ -120960.0_real64, 1814400.0_real64, -13305600.0_real64, &
            66528000.0_real64, -299376000.0_real64, -148262400.0_real64, 558835200.0_real64, &
            -239500800.0_real64, 99792000.0_real64, -33264000.0_real64, 7983360.0_real64, &
            -1209600.0_real64, 86400.0_real64/
        data (derivative_numerators(7,weight_index),weight_index=1,13)/ 86400.0_real64, -1244160.0_real64, 8553600.0_real64, &
            -38016000.0_real64, 128304000.0_real64, -410572800.0_real64, 0.0_real64, &
            410572800.0_real64, -128304000.0_real64, 38016000.0_real64, -8553600.0_real64, &
            1244160.0_real64,  - 86400.0_real64/
        data (derivative_numerators(8,weight_index),weight_index=1,13)/ -86400.0_real64, 1209600.0_real64, -7983360.0_real64, &
            33264000.0_real64, -99792000.0_real64, 239500800.0_real64, -558835200.0_real64, &
            148262400.0_real64, 299376000.0_real64, -66528000.0_real64, 13305600.0_real64, &
            -1814400.0_real64, 120960.0_real64/
        data (derivative_numerators(9,weight_index),weight_index=1,13)/ 120960.0_real64, -1658880.0_real64, 10644480.0_real64, &
            -42577920.0_real64, 119750400.0_real64, -255467520.0_real64, 447068160.0_real64, &
            -766402560.0_real64, 303937920.0_real64, 212889600.0_real64, -31933440.0_real64, &
            3870720.0_real64,  - 241920.0_real64/
        data (derivative_numerators(10,weight_index),weight_index=1,13)/ -241920.0_real64, 3265920.0_real64, -20528640.0_real64, &
            79833600.0_real64, -215550720.0_real64, 431101440.0_real64, -670602240.0_real64, &
            862202880.0_real64, -1077753600.0_real64, 476910720.0_real64, 143700480.0_real64, &
            -13063680.0_real64, 725760.0_real64/
        data (derivative_numerators(11,weight_index),weight_index=1,13)/ 725760.0_real64, -9676800.0_real64, 59875200.0_real64, &
            -228096000.0_real64, 598752000.0_real64, -1149603840.0_real64, 1676505600.0_real64, &
            -1916006400.0_real64, 1796256000.0_real64, -1596672000.0_real64, 684478080.0_real64, &
            87091200.0_real64,  - 3628800.0_real64/
        data (derivative_numerators(12,weight_index),weight_index=1,13)/ -3628800.0_real64, 47900160.0_real64, &
            & -292723200.0_real64, &
            1097712000.0_real64, -2822688000.0_real64, 5269017600.0_real64, -7376624640.0_real64, &
            7903526400.0_real64, -6586272000.0_real64, 4390848000.0_real64, -2634508800.0_real64, &
            967524480.0_real64, 39916800.0_real64/
        data (derivative_numerators(13,weight_index),weight_index=1,13)/ 39916800.0_real64, -522547200.0_real64, &
            & 3161410560.0_real64&
            , -11708928000.0_real64, 29638224000.0_real64, -54195609600.0_real64, &
            73766246400.0_real64, -75873853440.0_real64, 59276448000.0_real64, &
            -35126784000.0_real64, 15807052800.0_real64, -5748019200.0_real64, 1486442880.0_real64&
            /
        data derivative_denominator/ 479001600.0_real64/
        data starting_weights/ 1771.0_real64, 9235.0_real64, 5890.0_real64, 4610.0_real64, 35.0_real64, 59.0_real64/
        data starting_denominator/ 5760.0_real64/
        data (five_point_numerators(weight_index,2),weight_index=1,5)/ 251.0_real64, 646.0_real64, -264.0_real64, &
            & 106.0_real64,  - &
            19.0_real64/
        data (five_point_numerators(weight_index,3),weight_index=1,5)/ 232.0_real64, 992.0_real64, 192.0_real64, &
            & 32.0_real64,  - 8.0_real64&
            /
        data (five_point_numerators(weight_index,4),weight_index=1,5)/ 243.0_real64, 918.0_real64, 648.0_real64, &
            & 378.0_real64,  - &
            27.0_real64/
        data (five_point_numerators(weight_index,5),weight_index=1,5)/ 224.0_real64, 1024.0_real64, 384.0_real64, 1024.0_real64, &
            224.0_real64/
        data five_point_denominator/ 720.0_real64/
        data (six_point_numerators(weight_index,2),weight_index=1,6)/ 475.0_real64, 1427.0_real64, -798.0_real64, 482.0_real64, &
            -173.0_real64, 27.0_real64/
        data (six_point_numerators(weight_index,3),weight_index=1,6)/ 448.0_real64, 2064.0_real64, 224.0_real64, 224.0_real64, &
            -96.0_real64, 16.0_real64/
        data (six_point_numerators(weight_index,4),weight_index=1,6)/ 459.0_real64, 1971.0_real64, 1026.0_real64, 1026.0_real64, &
            -189.0_real64, 27.0_real64/
        data (six_point_numerators(weight_index,5),weight_index=1,6)/ 448.0_real64, 2048.0_real64, 768.0_real64, 2048.0_real64, &
            448.0_real64, 0.0_real64/
        data (six_point_numerators(weight_index,6),weight_index=1,6)/ 475.0_real64, 1875.0_real64, 1250.0_real64, 1250.0_real64, &
            1875.0_real64, 475.0_real64/
        data six_point_denominator/ 1440.0_real64/
        data initialize_derivative_table/ .true./
        ! The unit-spaced derivative weights are independent of h. The integration weights below scale with h.
        if (initialize_derivative_table) then
            weight_scale = 1.0_real64/derivative_denominator
            do column_index = 1,13
                do weight_index = 1,13
                    derivative_weights(weight_index,column_index) = derivative_numerators(weight_index,column_index)*weight_scale
                end do
            end do
            initialize_derivative_table = .false.
        end if
        starting_step_weights(1) = starting_weights(1)/starting_denominator
        weight_scale = transformed_grid_step/starting_denominator
        do weight_index = 2,6
            starting_step_weights(weight_index) = starting_weights(weight_index)*weight_scale
        end do
        weight_scale = transformed_grid_step/five_point_denominator
        do column_index = 2,4
            do weight_index = 2,5
                five_point_weights(weight_index,column_index) = weight_scale*five_point_numerators(weight_index,column_index)
            end do
        end do
        endpoint_weight = weight_scale*five_point_numerators(1,5)
        odd_interior_weight = weight_scale*five_point_numerators(2,5)
        middle_weight = weight_scale*five_point_numerators(3,5)
        shared_endpoint_weight = endpoint_weight+endpoint_weight
        weight_scale = transformed_grid_step/six_point_denominator
        do column_index = 2,6
            do weight_index = 1,6
                six_point_weights(weight_index,column_index) = weight_scale*six_point_numerators(weight_index,column_index)
            end do
        end do
        return
    end subroutine setqic
end module radial_stencil_weights
