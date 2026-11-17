! Merge validated left/reference and right/current bases without text round trips.
module csf_basis_merge
   use iso_fortran_env, only: int8, int64
   use parameter_def, only: max_orbitals
   use def_c, only: electron_count, nelecr
   use orb_c, only: csf_count, ncfr, nwr, npr, nakr, nhr, iqar, orbital_count, &
                    orbital_principal_n, orbital_kappa, orbital_l, orbital_two_j, &
                    orbital_symmetry_labels, csf_subshell_occupations
   use stat_c, only: jqsar, jcupar, csf_subshell_quantum_codes, csf_coupling_quantum_codes
   implicit none (type, external)
   private
   public :: merge_csf_bases
contains
   subroutine merge_csf_bases(left_core,right_core,status,message)
      integer, intent(in) :: left_core,right_core
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      integer :: merged_n(max_orbitals), merged_k(max_orbitals)
      character(len=2) :: merged_labels(max_orbitals)
      integer :: left, right, choice, count, i, j, rows, allocation_status
      integer(int64) :: wide_rows
      logical :: left_common, right_common
      integer(int8), pointer :: occupations(:,:), quantum(:,:,:), couplings(:,:)
      nullify(occupations,quantum,couplings)
      status=1
      message="left and right electron counts or cores differ"
      if (electron_count /= nelecr .or. left_core /= right_core) return
      if (left_core < 0 .or. left_core > min(nwr,orbital_count)) return
      if (any(npr(:left_core) /= orbital_principal_n(:right_core)) .or. &
          any(nakr(:left_core) /= orbital_kappa(:right_core))) return
      message="input basis arrays are not loaded"
      if (.not.associated(iqar) .or. .not.associated(jqsar) .or. .not.associated(jcupar)) return
      if (.not.associated(csf_subshell_occupations) .or. .not.associated(csf_subshell_quantum_codes) .or. &
          .not.associated(csf_coupling_quantum_codes)) return

      ! Stable topological merge. Common orbitals must have the same order;
      ! unconstrained left-only orbitals retain priority over right-only ones.
      left=1
      right=1
      count=0
      do while (left <= nwr .or. right <= orbital_count)
         if (left > nwr) then
            choice=2
         else if (right > orbital_count) then
            choice=1
         else if (npr(left) == orbital_principal_n(right) .and. nakr(left) == orbital_kappa(right)) then
            choice=3
         else
            left_common=any(npr(left) == orbital_principal_n(right:orbital_count) .and. &
                            nakr(left) == orbital_kappa(right:orbital_count))
            right_common=any(orbital_principal_n(right) == npr(left:nwr) .and. &
                             orbital_kappa(right) == nakr(left:nwr))
            if (.not.left_common) then
               choice=1
            else if (.not.right_common) then
               choice=2
            else
               message="inconsistent orbital order between the two bases"
               return
            end if
         end if
         message="merged orbital count exceeds configured capacity"
         if (count == max_orbitals) return
         count=count+1
         if (choice == 2) then
            merged_n(count)=orbital_principal_n(right)
            merged_k(count)=orbital_kappa(right)
            merged_labels(count)=orbital_symmetry_labels(right)
         else
            merged_n(count)=npr(left)
            merged_k(count)=nakr(left)
            merged_labels(count)=nhr(left)
         end if
         if (choice /= 2) left=left+1
         if (choice /= 1) right=right+1
      end do
      wide_rows=int(ncfr,int64)+int(csf_count,int64)
      message="merged CSF count exceeds suite integer capacity"
      if (wide_rows > int(huge(rows),int64) .or. wide_rows < 1) return
      rows=int(wide_rows)
      message="cannot allocate merged CSF arrays"
      allocate(occupations(max_orbitals,rows),stat=allocation_status)
      if (allocation_status /= 0) goto 900
      allocate(quantum(max_orbitals,3,rows),stat=allocation_status)
      if (allocation_status /= 0) goto 900
      allocate(couplings(max_orbitals,rows),stat=allocation_status)
      if (allocation_status /= 0) goto 900
      occupations=0_int8
      quantum=0_int8
      quantum(:count,3,:)=1_int8
      couplings(:,:ncfr)=int(jcupar(:,:ncfr),int8)
      couplings(:,ncfr+1:)=csf_coupling_quantum_codes
      do i=1,count
         do j=1,nwr
            if (merged_n(i) /= npr(j) .or. merged_k(i) /= nakr(j)) cycle
            occupations(i,:ncfr)=int(iqar(j,:ncfr),int8)
            quantum(i,:,:ncfr)=int(jqsar(j,:,:ncfr),int8)
         end do
         do j=1,orbital_count
            if (merged_n(i) /= orbital_principal_n(j) .or. merged_k(i) /= orbital_kappa(j)) cycle
            occupations(i,ncfr+1:)=csf_subshell_occupations(j,:)
            quantum(i,:,ncfr+1:)=csf_subshell_quantum_codes(j,:,:)
         end do
      end do
      deallocate(csf_subshell_occupations,csf_subshell_quantum_codes,csf_coupling_quantum_codes)
      csf_subshell_occupations=>occupations
      csf_subshell_quantum_codes=>quantum
      csf_coupling_quantum_codes=>couplings
      orbital_count=count
      csf_count=rows
      orbital_principal_n(:count)=merged_n(:count)
      orbital_kappa(:count)=merged_k(:count)
      orbital_symmetry_labels(:count)=merged_labels(:count)
      orbital_two_j(:count)=2*abs(merged_k(:count))-1
      where (merged_k(:count) > 0)
         orbital_l(:count)=merged_k(:count)
      elsewhere
         orbital_l(:count)=-merged_k(:count)-1
      end where
      status=0
      message=""
      return
900   continue
      if (associated(occupations)) deallocate(occupations)
      if (associated(quantum)) deallocate(quantum)
      if (associated(couplings)) deallocate(couplings)
   end subroutine merge_csf_bases
end module csf_basis_merge
