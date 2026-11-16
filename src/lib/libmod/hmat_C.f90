      module hmat_c
      use, intrinsic :: iso_fortran_env, only: real64, int64, double => real64, long => int64
      implicit none (type, external)
      private :: real64, int64
!...Created by Pacific-Sierra Research 77to90  4.3E  06:33:54  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer(int64) :: stored_matrix_element_count
      real(real64), pointer :: hamiltonian_values(:)
      integer, pointer :: stored_segment_end_offsets(:)
      integer, pointer :: matrix_partner_indices(:)
      integer :: ntpitmp(6)
      integer(int64) :: nelmnttmp
      integer :: ncoeitmp, ncoectmp, ncteitmp, nctectmp, nmcbptmp,  &
                 ncoretmp, nvpitmp, nkeitmp, nvintitmp, ncftmp
      real(real64) :: cutofftmp
      end module hmat_c
