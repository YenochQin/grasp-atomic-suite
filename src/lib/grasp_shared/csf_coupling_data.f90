! Stores subshell quantum numbers and cumulative coupling codes for the active and paired CSF bases.
      module csf_coupling_data
      use, intrinsic :: iso_fortran_env, only: int8
      implicit none (type, external)
      private :: int8
      integer, pointer :: jqsar(:, :, :) => null()
      integer(int8), pointer :: csf_subshell_quantum_codes(:, :, :) => null()
      integer, pointer :: jcupar(:, :) => null()
      integer(int8), pointer :: csf_coupling_quantum_codes(:, :) => null()
      end module csf_coupling_data
