module davidson_orthonormalization
    implicit none (type, external)
    private
    public :: mgs_nrm
contains
      subroutine mgs_nrm(vector_length, existing_vector_count, new_vector_count, projection_coefficients, &
         vectors)
      ! Apply two passes of modified Gram-Schmidt to the appended columns.
      use, intrinsic :: iso_fortran_env, only: real64
      use blas_lapack_interfaces, only: ddot, dgemv, dscal
      implicit none (type, external)
      integer, intent(in) :: vector_length
      integer , intent(in) :: existing_vector_count
      integer, intent(in) :: new_vector_count
      real(real64), intent(inout) :: projection_coefficients(new_vector_count)
      real(real64), intent(inout) :: vectors((existing_vector_count + new_vector_count)*vector_length)
      integer :: pass_index, new_vectors_start, reference_vector_start, reference_vector_index, &
         target_vector_start, target_vector_index, component_offset
      real(real64) :: vector_norm

      do pass_index = 1, 2

         new_vectors_start = existing_vector_count*vector_length + 1
         reference_vector_start = 1
         ! Remove projections onto the existing constraints and basis columns.
         do reference_vector_index = 1, existing_vector_count
            target_vector_start = new_vectors_start
            call dgemv ("T", vector_length, new_vector_count, 1.0_real64, vectors(target_vector_start), &
               vector_length, vectors(reference_vector_start), 1, 0.0_real64, projection_coefficients, 1)
            do target_vector_index = 1, new_vector_count
               do component_offset = 0, vector_length - 1
                  vectors(target_vector_start+component_offset) &
                     = vectors(target_vector_start+component_offset) &
                     - projection_coefficients(target_vector_index) &
                     *vectors(reference_vector_start+component_offset)
               end do
               target_vector_start = target_vector_start + vector_length
            end do
            reference_vector_start = reference_vector_start + vector_length
         end do
         ! Normalize each new column before projecting the subsequent new columns against it.
         do reference_vector_index = 1, new_vector_count
            target_vector_start = reference_vector_start + vector_length
            vector_norm = ddot(vector_length,vectors(reference_vector_start),1,vectors(reference_vector_start),1)
            vector_norm = sqrt(vector_norm)
            call dscal (vector_length, 1/vector_norm, vectors(reference_vector_start), 1)
            if (reference_vector_index < new_vector_count) then
            call dgemv ("T", vector_length, new_vector_count - reference_vector_index, 1.0_real64, &
               vectors(target_vector_start), vector_length, vectors(reference_vector_start), 1, 0.0_real64, &
               projection_coefficients, 1)
            end if
            do target_vector_index = reference_vector_index + 1, new_vector_count
               do component_offset = 0, vector_length - 1
                  vectors(target_vector_start+component_offset) &
                     = vectors(target_vector_start+component_offset) &
                     - projection_coefficients(target_vector_index-reference_vector_index) &
                     *vectors(reference_vector_start+component_offset)
               end do
               target_vector_start = target_vector_start + vector_length
            end do
            reference_vector_start = reference_vector_start + vector_length
         end do
      end do

      return
      end subroutine mgs_nrm
end module davidson_orthonormalization
