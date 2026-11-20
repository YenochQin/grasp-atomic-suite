! Dispatch a pair transfer that changes four distinct occupations.
module two_electron_four_shell_transfer
    implicit none (type, external)
    private
    public :: el5
contains
      subroutine el5(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,interaction_type)
      use csf_pair_angular_state,          only: peel_orbital_count => npeel
      use four_shell_transfer_12_34, only: el51
      use four_shell_transfer_13_24, only: el52
      use four_shell_transfer_14_23, only: el53
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
          interaction_type
      if(peel_orbital_count <= 3)return
      if(bra_shell_2 < ket_shell_1) then
        call el51(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,1,interaction_type)
      else if(bra_shell_1 > ket_shell_2.and.bra_shell_2 > ket_shell_2) then
        call el51(bra_csf_index,ket_csf_index,ket_shell_1,ket_shell_2,bra_shell_1,bra_shell_2,2,interaction_type)
      else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 < ket_shell_2.and.bra_shell_1 < ket_shell_1) then
        call el52(bra_csf_index,ket_csf_index,bra_shell_1,ket_shell_1,bra_shell_2,ket_shell_2,1,interaction_type)
      else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 > ket_shell_2.and.bra_shell_1 > ket_shell_1) then
        call el52(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,ket_shell_2,bra_shell_2,2,interaction_type)
      else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 > ket_shell_2.and.bra_shell_1 < ket_shell_1) then
        call el53(bra_csf_index,ket_csf_index,bra_shell_1,ket_shell_1,ket_shell_2,bra_shell_2,1,interaction_type)
      else if(bra_shell_2 > ket_shell_1.and.bra_shell_2 < ket_shell_2.and.bra_shell_1 > ket_shell_1) then
        call el53(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,bra_shell_2,ket_shell_2,2,interaction_type)
      else
        write(99,'(5x,"ERRO IN EL5 ")')
        stop
      end if
      return
      end subroutine el5
end module two_electron_four_shell_transfer
