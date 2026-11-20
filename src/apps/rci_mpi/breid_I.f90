      module breid_i
      implicit none (type, external)
      interface
      subroutine breid (bra_csf_index, ket_csf_index, ja1, ipca, jb1)
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index
      integer, intent(in) :: ket_csf_index
      integer, intent(in) :: ja1
      integer, intent(in) :: ipca
      integer, intent(in) :: jb1
      end subroutine breid
      end interface
      end module breid_i
