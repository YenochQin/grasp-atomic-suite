! Stores the number and values of parameters for the nuclear charge model.
      module nuclear_model_parameters_state
      use, intrinsic :: iso_fortran_env, only: real64
      implicit none (type, external)
      private :: real64
      integer :: nuclear_model_parameter_count
      real(real64) :: nuclear_model_parameters(2)
      end module nuclear_model_parameters_state
