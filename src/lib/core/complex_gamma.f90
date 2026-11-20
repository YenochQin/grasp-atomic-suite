! Evaluate the real and imaginary parts of Gamma using recurrence and a Stirling expansion.
module complex_gamma
    implicit none (type, external)
    private
    public :: cgamma
contains
      subroutine cgamma(argument_real, argument_imaginary, gamma_real, gamma_imaginary)
      use, intrinsic :: iso_fortran_env, only: real64
      use atomic_calculation_state, only: &
          maximum_exp_argument => expmax, minimum_exp_argument => expmin, pi, machine_precision => precis
      use nonnegative_atan2, only: arctan
      implicit none (type, external)
      real(real64), intent(in) :: argument_imaginary, argument_real
      real(real64), intent(out) :: gamma_real, gamma_imaginary
      integer :: stirling_term_index
      real(real64) :: stirling_coefficients(7), bernoulli_denominators(7)
      real(real64) :: half_log_two_pi, twice_term_index, distance_to_integer, argument_magnitude, log_gamma_imaginary, &
         & log_reflection_factor, shifted_argument_real, recurrence_product, &
         log_gamma_real, scale_factor, inverse_argument_squared, shifted_argument_imaginary, imaginary_part_squared, &
            & recurrence_product_real, recurrence_product_imaginary, term_real, term_imaginary, real_part_squared, &
         inverse_square_real, inverse_square_imaginary, argument_power_real, argument_power_imaginary
      logical,save :: initialize_stirling_coefficients, use_reflection_formula
      data stirling_coefficients/ 1.0_real64, -1.0_real64, 1.0_real64, -1.0_real64, 5.0_real64, -691.0_real64, 7.0_real64/
      data bernoulli_denominators/ 6.0_real64, 30.0_real64, 42.0_real64, 30.0_real64, 66.0_real64, 2730.0_real64, 6.0_real64/
      data half_log_two_pi/ 1.0_real64/
      data initialize_stirling_coefficients/ .true./
      ! Compute the seven scaled Bernoulli coefficients once for both real and complex arguments.
      if (initialize_stirling_coefficients) then
         half_log_two_pi = 0.5_real64*log(pi + pi)
         do stirling_term_index = 1, 7
            stirling_coefficients(stirling_term_index) = &
               & stirling_coefficients(stirling_term_index)/bernoulli_denominators(stirling_term_index)
            twice_term_index = dble(stirling_term_index + stirling_term_index)
            stirling_coefficients(stirling_term_index) = &
               & stirling_coefficients(stirling_term_index)/(twice_term_index*(twice_term_index - 1.0_real64))
         end do
         initialize_stirling_coefficients = .false.
      end if
      if (argument_imaginary == 0.0_real64) then
         ! Reflect negative real arguments and reject values close to an integer pole.
         if (argument_real <= 0.0_real64) then
            distance_to_integer = abs(dble(nint(argument_real)) - argument_real)
            if (distance_to_integer <= machine_precision + machine_precision) then
               write (*, '("CGAMMA: Argument (",1p,1d19.12,",",1d19.12,")"," too close to a pole.")') &
                  argument_real, argument_imaginary
               stop
            else
               argument_magnitude = pi/(-argument_real*sin(pi*argument_real))
               if (argument_magnitude < 0.0_real64) then
                  argument_magnitude = -argument_magnitude
                  log_gamma_imaginary = pi
               else
                  log_gamma_imaginary = 0.0_real64
               end if
               log_reflection_factor = log(argument_magnitude)
               shifted_argument_real = -argument_real
               use_reflection_formula = .true.
            end if
         else
            log_gamma_imaginary = 0.0_real64
            shifted_argument_real = argument_real
            use_reflection_formula = .false.
         end if
         recurrence_product = 1.0_real64
    2    continue
         if (shifted_argument_real < 10.0_real64) then
            recurrence_product = recurrence_product*shifted_argument_real
            shifted_argument_real = shifted_argument_real + 1.0_real64
            goto 2
         end if
         log_gamma_real = (shifted_argument_real - 0.5_real64)*log(shifted_argument_real) - shifted_argument_real + &
            & half_log_two_pi
         scale_factor = shifted_argument_real
         inverse_argument_squared = 1.0_real64/(shifted_argument_real*shifted_argument_real)
         do stirling_term_index = 1, 7
            scale_factor = scale_factor*inverse_argument_squared
            log_gamma_real = log_gamma_real + stirling_coefficients(stirling_term_index)*scale_factor
         end do
         log_gamma_real = log_gamma_real - log(recurrence_product)
         if (use_reflection_formula) log_gamma_real = log_reflection_factor - log_gamma_real
      else
         shifted_argument_real = argument_real
         shifted_argument_imaginary = argument_imaginary
         imaginary_part_squared = shifted_argument_imaginary*shifted_argument_imaginary
         recurrence_product_real = 1.0_real64
         recurrence_product_imaginary = 0.0_real64
    ! Shift the complex argument until its magnitude is large enough for the log-Gamma series.
    4    continue
         argument_magnitude = sqrt(shifted_argument_real*shifted_argument_real + imaginary_part_squared)
         if (argument_magnitude < 10.0_real64) then
            term_real = recurrence_product_real*shifted_argument_real - recurrence_product_imaginary*shifted_argument_imaginary
            term_imaginary = recurrence_product_real*shifted_argument_imaginary + &
               & recurrence_product_imaginary*shifted_argument_real
            recurrence_product_real = term_real
            recurrence_product_imaginary = term_imaginary
            shifted_argument_real = shifted_argument_real + 1.0_real64
            goto 4
         end if
         real_part_squared = shifted_argument_real*shifted_argument_real
         term_real = 0.5_real64*log(real_part_squared + imaginary_part_squared)
         term_imaginary = arctan(shifted_argument_imaginary,shifted_argument_real)
         log_gamma_real = (shifted_argument_real - 0.5_real64)*term_real - shifted_argument_imaginary*term_imaginary - &
            & shifted_argument_real + half_log_two_pi
         log_gamma_imaginary = (shifted_argument_real - 0.5_real64)*term_imaginary + shifted_argument_imaginary*term_real &
            & - shifted_argument_imaginary
         scale_factor = (real_part_squared + imaginary_part_squared)**(-2)
         inverse_square_real = (real_part_squared - imaginary_part_squared)*scale_factor
         inverse_square_imaginary = -2.0_real64*shifted_argument_real*shifted_argument_imaginary*scale_factor
         argument_power_real = shifted_argument_real
         argument_power_imaginary = shifted_argument_imaginary
         do stirling_term_index = 1, 7
            term_real = argument_power_real*inverse_square_real - argument_power_imaginary*inverse_square_imaginary
            term_imaginary = argument_power_real*inverse_square_imaginary + argument_power_imaginary*inverse_square_real
            scale_factor = stirling_coefficients(stirling_term_index)
            log_gamma_real = log_gamma_real + term_real*scale_factor
            log_gamma_imaginary = log_gamma_imaginary + term_imaginary*scale_factor
            argument_power_real = term_real
            argument_power_imaginary = term_imaginary
         end do
         log_gamma_real = log_gamma_real - 0.5_real64*log(recurrence_product_real*recurrence_product_real + &
            & recurrence_product_imaginary*recurrence_product_imaginary)
         log_gamma_imaginary = log_gamma_imaginary - arctan(recurrence_product_imaginary,recurrence_product_real)
      end if
      if (log_gamma_real<=maximum_exp_argument .and. log_gamma_real>=minimum_exp_argument) then
         scale_factor = exp(log_gamma_real)
      else
         write (*, '("CGAMMA: Argument to exponential function"," (",1p,1d19.12,") out of range.")') log_gamma_real
         stop
      end if
      gamma_real = scale_factor*cos(log_gamma_imaginary)
      gamma_imaginary = scale_factor*sin(log_gamma_imaginary)
      return
      return
      end subroutine cgamma
end module complex_gamma
