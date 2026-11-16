      module pos_c
      use vast_kind_param, only:  double
!...Created by Pacific-Sierra Research 77to90  4.3E  06:38:40  12/28/06
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
      integer :: total_csf_count, mixing_coefficient_storage_size
      integer, dimension(:), pointer :: csf_offsets_by_block
      integer, dimension(:), pointer :: eigenvalue_offsets_by_block, mixing_offsets_by_block

      end module pos_c
