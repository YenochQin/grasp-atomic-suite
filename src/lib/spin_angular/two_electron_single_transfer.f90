! Dispatch two-electron terms with a net transfer of one electron.
module two_electron_single_transfer
    implicit none (type, external)
    private
    public :: el3
contains
      subroutine el3(bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2,interaction_type)
      use csf_pair_angular_state,         only:  peel_orbital_count => npeel
      use single_transfer_triple_on_first_shell, only: el31
      use single_transfer_triple_on_second_shell, only: el32
      use single_transfer_with_spectator_shell, only: el33
      implicit none (type, external)
      integer, intent(in) :: bra_csf_index,ket_csf_index,bra_shell_1,bra_shell_2,ket_shell_1,ket_shell_2, &
          interaction_type
      if(peel_orbital_count <= 1)return
      if(bra_shell_2 == ket_shell_2) then
        if(bra_shell_1 == bra_shell_2.or.ket_shell_1 == bra_shell_2) then
          if(bra_shell_1 == ket_shell_1)goto 10
          if(ket_shell_1 /= bra_shell_2) then
            call el32(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,bra_shell_1,bra_shell_2,ket_shell_1, &
                ket_shell_2,interaction_type)
          else
            call el31(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,bra_shell_1,bra_shell_2,ket_shell_1, &
                ket_shell_2,interaction_type)
          end if
        else
          call el33(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_1,bra_shell_2,1,bra_shell_1,bra_shell_2, &
              ket_shell_1,ket_shell_2,interaction_type)
        end if
        return
      else if(bra_shell_1 == ket_shell_1) then
        if(bra_shell_2 == bra_shell_1.or.ket_shell_2 == bra_shell_1) then
          if(bra_shell_2 == ket_shell_2)goto 10
          if(ket_shell_2 /= bra_shell_1) then
            call el32(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_2,bra_shell_1,bra_shell_2,ket_shell_1, &
                ket_shell_2,interaction_type)
          else
            call el31(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_2,bra_shell_1,bra_shell_2,ket_shell_1, &
                ket_shell_2,interaction_type)
          end if
        else
          call el33(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_2,bra_shell_1,1,bra_shell_1,bra_shell_2, &
              ket_shell_1,ket_shell_2,interaction_type)
        end if
        return
      else if(bra_shell_1 == ket_shell_2) then
        if(bra_shell_2 == bra_shell_1.or.ket_shell_1 == bra_shell_1) then
          if(bra_shell_2 == ket_shell_1)goto 10
          if(ket_shell_1 /= ket_shell_2) then
             call el32(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_2,bra_shell_1,bra_shell_2,ket_shell_1, &
                 ket_shell_2,interaction_type)
          else
             call el31(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_2,bra_shell_1,bra_shell_2,ket_shell_1, &
                 ket_shell_2,interaction_type)
          end if
        else
          call el33(bra_csf_index,ket_csf_index,ket_shell_1,bra_shell_2,bra_shell_1,2,bra_shell_1,bra_shell_2, &
              ket_shell_2,ket_shell_1,interaction_type)
        end if
        return
      else if(bra_shell_2 == ket_shell_1) then
        if(bra_shell_1 == bra_shell_2.or.ket_shell_2 == bra_shell_2) then
          if(bra_shell_1 == ket_shell_2)goto 10
          if(ket_shell_2 /= bra_shell_2) then
            call el32(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_1,bra_shell_1,bra_shell_2,ket_shell_1, &
                ket_shell_2,interaction_type)
          else
            call el31(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_1,bra_shell_1,bra_shell_2,ket_shell_1, &
                ket_shell_2,interaction_type)
          end if
        else
          call el33(bra_csf_index,ket_csf_index,ket_shell_2,bra_shell_1,bra_shell_2,2,bra_shell_1,bra_shell_2, &
              ket_shell_2,ket_shell_1,interaction_type)
        end if
        return
      end if
   10 write(99,'(5x,"ERRO IN EL3  PMGG RAGG")')
      stop
      end subroutine el3
end module two_electron_single_transfer
