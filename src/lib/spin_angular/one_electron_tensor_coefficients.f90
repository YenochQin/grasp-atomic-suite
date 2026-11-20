! Generate one-electron angular coefficients for a specified rank and parity.
module one_electron_tensor_coefficients
    implicit none (type, external)
    private
    public :: oneparticlejj
contains
      subroutine oneparticlejj(tensor_rank,operator_parity,bra_csf_index,ket_csf_index,bra_orbital_index, &
          ket_orbital_index,shell_coefficients)
      use, intrinsic :: iso_fortran_env, only: real64
      use grasp_parameters,   only:  max_orbitals
      use numeric_constants, only: eleven, angular_tolerance => eps, four, half, one, seven, ten, tenth, three, two, zero
      use csf_pair_angular_state, only: bra_intermediate_couplings => jjc1, ket_intermediate_couplings => jjc2, &
          bra_subshell_quantum_codes => jjq1, ket_subshell_quantum_codes => jjq2, peel_orbital_indices => jlist, &
          core_orbital_indices => klist, max_orbitals, core_orbital_count => ncore, peel_orbital_count => npeel, &
          bra_occupations => nq1, ket_occupations => nq2
      use orbital_csf_basis,   only: orbital_count, orbital_kappa
      use angular_recoupling_backup,  only: saved_peel_orbitals => jlis, saved_bra_j_weights => jc1s, saved_ket_j_weights => jc2s
      use angular_momentum_triangle, only: itrig
      use one_electron_csf_filter, only: ichkq1
      use subshell_occupation_class, only: ichop
      use csf_quantum_queries, only: ispar, itjpo
      use csf_pair_quantum_numbers, only: setqna
      use one_electron_tensor_same_subshell, only: oneparticlejj1
      use one_electron_tensor_between_subshells, only: oneparticlejj2
      implicit none (type, external)
      integer, intent(in)  :: tensor_rank,operator_parity,bra_csf_index,ket_csf_index
      integer, intent(out) :: bra_orbital_index,ket_orbital_index
      real(real64), intent(out) :: shell_coefficients(max_orbitals)
      integer :: two_rank_plus_one,orbital_index,total_occupation_difference,bra_shell_position,ket_shell_position, &
          peel_position,occupation_difference,active_shell_count,diagonal_orbital,shell_position,shift_step, &
          first_orbital,shift_count,  &
                 two_bra_j_plus_one,two_ket_j_plus_one,orbital_difference,peel_coupling_count
      integer :: unused_orbital_labels(2),unused_orbital_j_weights(2)
      real(real64) :: angular_coefficient

      ! A zero bra orbital index signals that no coefficient channel has been identified.
      bra_orbital_index = 0
      two_rank_plus_one = tensor_rank+tensor_rank+1
      if(itrig(itjpo(bra_csf_index),itjpo(ket_csf_index),two_rank_plus_one)==0)return
      if((operator_parity/=0).and.(ispar(bra_csf_index)*ispar(ket_csf_index)*operator_parity/=1))return
      if(ichkq1(bra_csf_index,ket_csf_index)==0)return
      call setqna (bra_csf_index,ket_csf_index)
      do orbital_index = 1,orbital_count
         shell_coefficients(orbital_index) = zero
      end do
      total_occupation_difference = 0
      bra_shell_position = 0
      ket_shell_position = 0
      if (peel_orbital_count /= 0) then
        do peel_position = 1,peel_orbital_count
          orbital_index = peel_orbital_indices(peel_position)
          occupation_difference = bra_occupations(orbital_index)-ket_occupations(orbital_index)
          if (abs (occupation_difference) > 1) return
          if (occupation_difference > 0) then
            bra_shell_position = peel_position
            total_occupation_difference = total_occupation_difference+1
          else if (occupation_difference < 0) then
            ket_shell_position = peel_position
            total_occupation_difference = total_occupation_difference+1
          end if
        end do
        if (total_occupation_difference > 2) return
        active_shell_count = peel_orbital_count
      end if

      ! A one-electron transfer fixes both active subshells. Equal occupations require a scan.
      if (total_occupation_difference == 2) goto 19
      diagonal_orbital = 0
      if (peel_orbital_count == 0) goto 9
      do shell_position = 1,peel_orbital_count
         saved_peel_orbitals(shell_position) = peel_orbital_indices(shell_position)
      end do
      if (peel_orbital_count == 1) goto 9
      peel_coupling_count = peel_orbital_count-1
      do shell_position = 1,peel_coupling_count
         saved_bra_j_weights(shell_position) = bra_intermediate_couplings(shell_position)
         saved_ket_j_weights(shell_position) = ket_intermediate_couplings(shell_position)
      end do
    9 diagonal_orbital = diagonal_orbital+1
      if (diagonal_orbital > orbital_count) return
      if (ichop (diagonal_orbital,bra_csf_index) == -1) goto 9
      if (ichop (diagonal_orbital,bra_csf_index) == 0) goto 16
      shell_position = 1
      if (peel_orbital_count==0) goto 15
      do shell_position = 1,peel_orbital_count
        orbital_index = peel_orbital_indices(shell_position)
        if (diagonal_orbital < orbital_index) goto 11
      end do
      shell_position = peel_orbital_count+1
      goto 13

   ! Open a slot for a closed subshell in the ordered active coupling chain.
   11 shift_count = peel_orbital_count-shell_position+1
      do shift_step = 1,shift_count
         peel_orbital_indices(peel_orbital_count+2-shift_step) = peel_orbital_indices(peel_orbital_count+1-shift_step)
         if (peel_orbital_count==shift_step) goto 13
         bra_intermediate_couplings(peel_orbital_count+1-shift_step) = bra_intermediate_couplings( &
             peel_orbital_count-shift_step)
         ket_intermediate_couplings(peel_orbital_count+1-shift_step) = ket_intermediate_couplings( &
             peel_orbital_count-shift_step)
      end do
   13 continue
      if (shell_position < 3) goto 14
      bra_intermediate_couplings(shell_position-1) = bra_intermediate_couplings(shell_position-2)
      ket_intermediate_couplings(shell_position-1) = ket_intermediate_couplings(shell_position-2)
      goto 15
   14 first_orbital = peel_orbital_indices(1)
      bra_intermediate_couplings(1) = bra_subshell_quantum_codes(3,first_orbital)
      ket_intermediate_couplings(1) = ket_subshell_quantum_codes(3,first_orbital)
   15 peel_orbital_indices(shell_position) = diagonal_orbital
      bra_shell_position = shell_position
      ket_shell_position = shell_position
      active_shell_count = peel_orbital_count+1
      goto 19
   16 active_shell_count = peel_orbital_count
      do peel_position = 1,peel_orbital_count
        orbital_difference = diagonal_orbital-peel_orbital_indices(peel_position)
        if (orbital_difference==0) goto 18
      end do
   18 bra_shell_position = peel_position
      ket_shell_position = peel_position
   19 bra_orbital_index = peel_orbital_indices(bra_shell_position)
      ket_orbital_index = peel_orbital_indices(ket_shell_position)
      two_bra_j_plus_one = 2*abs (orbital_kappa(bra_orbital_index))
      two_ket_j_plus_one = 2*abs (orbital_kappa(ket_orbital_index))
      if (itrig (two_bra_j_plus_one,two_ket_j_plus_one,two_rank_plus_one)==1) goto 99
      if (total_occupation_difference == 2) return
      goto 100
   99 continue
      if (total_occupation_difference == 0) then
          call oneparticlejj1(active_shell_count,tensor_rank,bra_csf_index,ket_csf_index,bra_shell_position, &
              ket_shell_position,angular_coefficient)
      else if (total_occupation_difference == 2) then
        call oneparticlejj2(active_shell_count,tensor_rank,bra_shell_position,ket_shell_position,angular_coefficient)

        ! Store a transfer coefficient in slot 1 and return the orbital indices separately.
        shell_coefficients(1) = angular_coefficient
        return
      end if

      ! An equal-orbital contribution is indexed by its full orbital label.
      shell_coefficients(diagonal_orbital) = angular_coefficient
100   continue
      if (peel_orbital_count == 0) goto 9
      do shell_position = 1,peel_orbital_count
         peel_orbital_indices(shell_position) = saved_peel_orbitals(shell_position)
      end do
      if (peel_orbital_count == 1) goto 9
      peel_coupling_count = peel_orbital_count-1
      do  shell_position = 1,peel_coupling_count
         bra_intermediate_couplings(shell_position)  = saved_bra_j_weights(shell_position)
         ket_intermediate_couplings(shell_position)  = saved_ket_j_weights(shell_position)
      end do
      goto 9
      return
      end subroutine oneparticlejj
end module one_electron_tensor_coefficients
