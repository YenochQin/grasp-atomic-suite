!
!***********************************************************************
!                                                                      *
      module jj2lsj_code
!                                                                      *
!     This module contains the procedures which are used to performe   *
!     the jj-LSJ transformation.                                       *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     Vilnius                                  last update: Jan 2017   *
!     Modified by Chongyang Chen                                2023   *
!     Some minor fixes labelled by CYC-2023                            *
!     Modified by G. Gaigalas                                   2024   *
!     Modified by G. Gaigalas                                   2025   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use jj2lsj_profile, only: profile_emit, profile_file_size, profile_finish, profile_start, profile_time
      use jj2lsj_postprocess, only: ls_ranking, post_worker_count, rank_ls
      implicit none (type, external)
      private :: real64
!-----------------------------------------------
!   R o u t i n e s
!-----------------------------------------------
      public  :: asf2ls
!               Expands an atomic state function, which is represented
!               in a jj-coupling CSF basis into a basis of LS-coupling CSF.
      public  :: coeflsjj
!               Returns the value of the LS-jj transformation matrix
!               for a given set of quantum numbers.
      private :: coeflsjj2
!               Returns the value of the LS-jj transformation matrix
!               (l^2 LSJ| j_1 j_2 J).
      private :: coeflsjjs
!               Returns the value of the LS-jj transformation matrix
!               form an array LS_jj_*_*.
      public  :: dallocasfls
!               Dellocates the storage of asf_set_LS.
      public  :: getchls
!               A spectroscopic notation of shell in LS coupling is return.
      public  :: getchls_lbl
!               A spectroscopic notation of shell in LS coupling is return.
      public  :: getxj
!
      public  :: gettermls
!               This procedure return all allowed subshell terms
!               (l, w, Q, L, S) for given l^N which must be 0, 1, 2 or 3.
      public  :: inscreen
!               The input from the screen.
      public  :: inscreenlev
!               Attempts to interpret the serial level numbers from a
!               string.
      public  :: jj2lsj
!               Controls the transformation of atomic states from a jj-
!               into a LS-coupled CSF basis.
      public  :: packlscsf
!                Encording all CSF in LS - coupling.
      public  :: prcsfjj
!               Prints all information about the single CSF in jj-coupling
      public  :: prcsfls
!               Prints all information about the single CSF in LS-coupling
      public  :: prcsflsall
!               Prints all information about the single CSF scheme csf_set
!               in LS-coupling
      public  :: setls  ! Basis construction for the conversion interface.
!               This subroutine fills up the variable asf_set_LS%csf_set_LS
!               with data generated using the one from asf_set%csf_set
!               .....................................................
!               This subroutine contains the following internal routines:
!                * subroutine setLS_action
!                  The subroutine defines the "action" of subroutine
!                  setLS_job_count: whether it counts
!                  the number of csfs_LS (asf_set_LS%csf_set_LS%novcsf)
!                  or fills the arrays of wave functions in LS coupling
!                  with asf_set_LS%csf_set_LS%csf(...) with
!                  the corresponding quantum nubers.
!                * subroutine setLS_add_quantum_numbers
!                  The subroutine adds quantum numbers stored
!                  in temprorary arrays Li, Si, L_i, S_i, w, Q to
!                  the corresponding arrays of asf_set_LS%csf_set_LS%csf().
!                  private :: setLS_job_count
!                * recursive subroutine setLS_job_count
!                  Recursive subroutine for the calculation of the
!                  number of csfs_LS and corresponding quantum numbers.
!                * function setLS_equivalent_csfs
!                  This subroutine defines the "equivalency" of two csfs_jj
!                  in the sence of generation of csfs_LS
!                  number of csfs_LS and corresponding quantum numbers.
!               .....................................................
      public  :: tralsjj
!               Return the value of the transformation matrix
!               from jj- to LS-coupling scheme in the case of any
!               number of open shells.
      public  :: tralsjjmp
!               Return the value of main part of the transformation
!               matrix from jj- to LS-coupling scheme in the
!               case of any number of open shells.
      private :: uniquelsj
!               Subroutine defines a unique labels for energy levels
!-----------------------------------------------
!   D e f i n i t i o n   o f   A S F
!   i n   L S -   C o u p l i n g
!-----------------------------------------------
      type, public :: nl
         integer :: n, l
      end type nl
!
      type, public :: cs_function_ls
         integer :: totalj
         character(len=1)  :: parity
         integer, pointer :: occupation(:) => null()
         integer, pointer :: seniority(:) => null()
         integer, pointer :: w(:) => null()
         integer, pointer :: shelll(:) => null()
         integer, pointer :: shells(:) => null()
         integer, pointer :: shelllx(:) => null()
         integer, pointer :: shellsx(:) => null()
      end type cs_function_ls
!
      type, public :: csf_basis_ls
         integer :: nocsf         ! Number of CSF in the basis.
         integer :: nwshells      ! Number of (nonrelativistic) shells.
         integer :: nwcore        ! Number of (closed) core shells.
         integer :: number_of_electrons
         type(nl), pointer             :: shell(:) => null()
         type(cs_function_ls), pointer :: csf(:) => null()
         type(parent_from_jj), pointer :: parent(:) => null()
      end type csf_basis_ls
!
      type, public :: as_function_ls
         integer          :: level_no
         integer          :: max_csf_no
         integer          :: totall, totals, totalj
         character(len=1) :: parity
         real(real64)     :: energy
         real(real64), pointer :: eigenvector(:) => null()
      end type as_function_ls
!
      type, public :: asf_basis_ls
         integer      :: noasf           ! Number of considered ASF.
         real(real64) :: average_energy  ! Averaged energy of this set of ASF.
         type(as_function_ls), pointer :: asf(:) => null()
         type(csf_basis_ls) :: csf_set_ls
      end type asf_basis_ls
!-----------------------------------------------
!   G l o b a l   V a r i a b l e s
!-----------------------------------------------
      type, public :: parent_from_jj
         integer :: parent_minus
         integer :: parent_plus
      end type parent_from_jj
!
      type, public :: lsj_list
         integer::list_size                    !number items in a list
         integer,pointer:: items(:) => null() !serial numbers of lists items
      end type lsj_list
!
      type(asf_basis_ls), public :: asf_set_ls
!
      integer, private :: imincompoff, ls_convention
      integer,pointer:: jcoup(:, :) ! xJ values
      real(real64) :: epsnew, mincomp
      character(len=1), private, parameter :: l_string(0:20) =&
               [ "S","P","D","F","G","H","I","K","L","M","N","O","Q",   &
                  "R","T","U","V","W","X","Y","Z" ]
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer, parameter :: blocks_number  = 20
      integer, parameter :: vectors_number = 100000
!-----------------------------------------------
!
contains
!
!***********************************************************************
!                                                                      *
!     Parallel rows are committed in the original LS order. The unchanged
!     scalar backend remains the reference and handles sequential *.lsj.T I/O.
      subroutine asf2ls(iw1,truncation_threshold_index,levmax,symmetry_block_index,levels,first_csf_index,last_csf_index,ioutt)
      use, intrinsic :: iso_fortran_env, only: real32
      use selected_level_indices, only: ivec
!$    USE omp_lib, ONLY: omp_get_max_threads
      implicit none (type, external)
      integer, intent(in) :: iw1, levmax, symmetry_block_index, first_csf_index, last_csf_index, ioutt
      integer, intent(in) :: truncation_threshold_index(:), levels(blocks_number,vectors_number)
      integer :: threads, batch_size, first, last, row, ls, lev, level, status
      integer, allocatable :: selected(:)
      real(real64), allocatable :: rows(:,:), cumulative(:)
      logical, allocatable :: active(:)
      logical :: finished
      character(len=16) :: serial_backend

      threads = 1
!$    threads = omp_get_max_threads()
      call get_environment_variable("JJ2LSJ_SERIAL", serial_backend, status=status)
      if (status == 0) then
         if (trim(serial_backend) == "1") threads = 1
      end if
      if (threads <= 1 .or. ioutt /= 0) then
         call asf2ls_serial(iw1,truncation_threshold_index,levmax,symmetry_block_index,levels,first_csf_index,last_csf_index,ioutt)
         return
      end if
      if (levmax <= 0 .or. asf_set_ls%csf_set_ls%nocsf <= 0) return

!     Bound scratch to about 8 MiB, even for many selected levels. Each row
!     is contiguous, with no per-pair allocations or cross-thread reduction.
      batch_size = min(256, 4*min(threads,64))
      batch_size = min(batch_size, max(1,1048576/levmax))
      batch_size = min(batch_size, asf_set_ls%csf_set_ls%nocsf)
      allocate(rows(levmax,batch_size), active(batch_size))
      allocate(cumulative(levmax), selected(levmax))
      selected = levels(symmetry_block_index,1:levmax)
      cumulative = 0.0_real64
      finished = .false.
      write(*,"(A,I0,A,I0)") " jj2lsj OpenMP: max threads=",threads, &
                           ", LS rows per batch=",batch_size

      do first = 1, asf_set_ls%csf_set_ls%nocsf, batch_size
         last = min(first+batch_size-1,asf_set_ls%csf_set_ls%nocsf)
!$omp parallel do default(none) schedule(dynamic,1) &
!$omp shared(first,last,iw1,truncation_threshold_index,selected,first_csf_index,last_csf_index,rows,active) private(row,ls)
         do ls = first, last
            row = ls-first+1
            call transform_ls_row(ls,iw1,truncation_threshold_index,selected,first_csf_index,last_csf_index, &
                                  rows(:,row),active(row))
         end do
!$omp end parallel do

!        Reproduce the scalar accumulation, strict comparison and first
!        stopping row exactly; speculative later rows are never committed.
         do ls = first, last
            if (mod(ls,2000) == 0) print *, "LS_number=",ls
            row = ls-first+1
            if (.not. active(row)) cycle
            do lev = 1,levmax
               level = selected(lev)
               cumulative(lev) = rows(lev,row)*rows(lev,row)+cumulative(lev)
               asf_set_ls%asf(level)%eigenvector(ls) = rows(lev,row)
               asf_set_ls%asf(level)%level_no = ivec(level)
            end do
            if (imincompoff == 1) then
               finished = .not. any(cumulative*100.0_real32+dabs(mincomp) <= 100.0_real32)
               if (finished) exit
            end if
         end do
         if (finished) exit
      end do
      deallocate(rows,active,cumulative,selected)
      end subroutine asf2ls


      recursive subroutine transform_ls_row(ls,iw1,ithresh,selected, &
                                           ncfmin,ncfmax,values,active)
      use eigenstate_expansions, only: csf_mixing_coefficients
      use orbital_csf_basis, only: csf_count
      use csf_quantum_queries, only: itjpo, ispar
      use, intrinsic :: iso_fortran_env, only: int64
      implicit none (type, external)
      integer, intent(in) :: ls, iw1, ncfmin, ncfmax
      integer, intent(in) :: ithresh(:), selected(:)
      real(real64), intent(out) :: values(:)
      logical, intent(out) :: active
      integer :: jj, lev, parity
      integer(int64) :: offset, ncf_long
      real(real64) :: coefficient

      values = 0.0_real64
      parity = ispar(ncfmin+iw1-1)
      active = (asf_set_ls%csf_set_ls%csf(ls)%parity == "+" .and. parity == 1) .or. &
               (asf_set_ls%csf_set_ls%csf(ls)%parity == "-" .and. parity == -1)
      if (.not. active) return
      ncf_long = csf_count
      do jj = ncfmin, ncfmax
         if (ithresh(jj) /= 1) cycle
         if (asf_set_ls%csf_set_ls%csf(ls)%totalj /= itjpo(jj)-1) cycle
         coefficient = tralsjj(jj,ls)
         do lev = 1,size(selected)
            offset = (selected(lev)-1)*ncf_long
            values(lev) = values(lev)+csf_mixing_coefficients(jj+offset)*coefficient
         end do
      end do
      end subroutine transform_ls_row


      subroutine asf2ls_serial(iw1,truncation_threshold_index,levmax,symmetry_block_index,levels,first_csf_index,     &
                        last_csf_index,ioutt)
      use, intrinsic :: iso_fortran_env, only: real32
!                                                                      *
!     Expands an atomic state functions from the same block,           *
!     which is represented in a jj-coupling CSF basis into a basis     *
!     of LS-coupling CSF.                                              *
!                                                                      *
!     Calls: ispar, itjpo, tranLSjj.                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!     Modified by C. Cychen                                     2021   *
!     Modified by G. Gaigalas                                   2022   *
!     Modified by G. Gaigalas                                   2024   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: int64
      use eigenstate_expansions,          only: csf_mixing_coefficients
      use selected_level_indices,          only: ivec
      use orbital_csf_basis,           only: csf_count
      use numeric_constants,          only: zero, one
!      USE JJ2LSJ_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: itjpo, ispar
!      USE idigit_I
!      USE ittk_I
!      USE ixjtik_I
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: iw1, levmax,symmetry_block_index,ioutt
      integer, intent(in) :: first_csf_index, last_csf_index
      integer, intent(in) :: truncation_threshold_index(:)
      integer, intent(in) :: levels(blocks_number, vectors_number)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=18) :: string
      integer      :: jj_number, lev, level, ls_number, ijj, ils
      integer      :: imincomp, idum
      real(real64) :: wa_transformation
      real(real64) :: wa(vectors_number)
      real(real64) :: wb(vectors_number)
      integer(int64) :: loc_i8, ncf_i8
!-----------------------------------------------
      wb = zero
      ncf_i8 = csf_count
      if(ioutt == 1) write(59) " *   Block Number=",symmetry_block_index
      if(ioutt == 2) then
         read(59) string, idum
         if(string(1:18) /= " *   Block Number=" .or.                  &
                                                  idum /= symmetry_block_index) then
            print*, "Error in transformation file *.lsj.T"
         end if
      end if
      do ls_number = 1, asf_set_ls%csf_set_ls%nocsf
         if(((ls_number/2000)*2000) == ls_number)                      &
                         print*, "LS_number=",ls_number
         if ((asf_set_ls%csf_set_ls%csf(ls_number)%parity == "+"       &
!CYC-2023 Corrections, perform transfermation for both even and odd
!states
            !.and.  ISPAR(iw1) == 1)  .or.                              &
            .and.  ispar(first_csf_index+iw1-1) == 1)  .or.                     &
            (asf_set_ls%csf_set_ls%csf(ls_number)%parity  == "-"       &
            !.and.  ISPAR(iw1) == -1)) then
            .and.  ispar(first_csf_index+iw1-1) == -1)) then
!CYC-2023 end
            wa = zero
            do jj_number = first_csf_index, last_csf_index
               if(truncation_threshold_index(jj_number) == 1 .and.                        &
                 (asf_set_ls%csf_set_ls%csf(ls_number)%totalj ==       &
                  itjpo(jj_number)-1)) then
                  if(ioutt <= 1)                                       &
                       wa_transformation = tralsjj(jj_number,ls_number)
                  if(ioutt == 1) write(59)                             &
                         jj_number-first_csf_index+1,ls_number,wa_transformation
                  if(ioutt == 2) then
                  read(59) ijj, ils, wa_transformation
                     if(ijj /= jj_number-first_csf_index+1 .or.                 &
                        ils /= ls_number ) then
                          print*, "Error in jj2lsj transformation file"
                          stop
                     end if
                  end if
                  do lev = 1,levmax
                     level = levels(symmetry_block_index,lev)
                     loc_i8 = (level-1)*ncf_i8
                     wa(lev) =                                         &
                      wa(lev)+csf_mixing_coefficients(jj_number+loc_i8)*wa_transformation
                  end do
               end if
            end do
            imincomp = 1
            do lev = 1,levmax
               level = levels(symmetry_block_index,lev)
               wb(lev) = wa(lev)*wa(lev) + wb(lev)
               asf_set_ls%asf(level)%eigenvector(ls_number) = wa(lev)
               asf_set_ls%asf(level)%level_no = ivec(level)
               if((wb(lev)*100.0_real32+dabs(mincomp)) <= 100.0_real32) imincomp = 0
            end do
            if(imincompoff == 1 ) then
               if(imincomp == 1) goto 1
            end if
         end if
      end do
    1 continue
      return
      end subroutine asf2ls_serial
!
!***********************************************************************
!                                                                      *
      recursive function coeflsjj(l_shell,n,w,q,l,s,j,jm_shell,nm,qm,jm,jp_shell,&
                                                     qp,jp) result(wa)
!                                                                      *
!     Returns the value of the LS-jj transformation matrix for a given *
!     set of quantum numbers form an array LS_jj_*_*.                  *
!                                                                      *
!     Note that all (generalized) angular momentum quantum numbers     *
!     except for l must be given as twice the original numbers, i.e.   *
!     for the quantum numbers Q, L, S, J, jm_shell, Qm, Jm, jp_shell,  *
!     Qp, Jp.                                                          *
!                                                                      *
!     Calls: coefLSjj2, coefLSjjs.                                     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use numeric_constants,          only: zero, one
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_shell, n, w, q, l, s, j,            &
                             jm_shell, nm, qm, jm, jp_shell, qp, jp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: wa
      integer      :: nn, n1, n2, factor
!-----------------------------------------------
      wa = zero; factor = one
      if (jm_shell < jp_shell .or. jm_shell == jp_shell)  then
         if (l_shell > 0 .and. n > 2*l_shell +1)  then
            nn = 4*l_shell + 2 - n;   n1 = jm_shell + 1 - nm
            n2 = jp_shell + 1 - n + nm
            if (mod(2*l_shell+1-q -((jm_shell+1)/2-qm)         &
                 -((jp_shell+1)/2-qp),4) /= 0) factor = - factor
         else
            nn = n;   n1 = nm;   n2 = n - nm
         end if
         if (nn == 1  .or.  nn == 0)  then
            if (nn == 0 .and. n1 == 0 .and. n2 == 0) then
               wa = one
            else if (n1 == 1 .and. n2 == 0 .and. j == jm) then
               wa = one
            else if (n1 == 0 .and. n2 == 1 .and. j == jp) then
               wa = one
            else
               wa = zero
            end if
         else if (nn == 2)  then
            if (j > jm + jp  .or.  j < abs(jm - jp)) then
               wa = zero
            else
               if (n1 == 2 .and. n2 == 0) then
                  wa = coeflsjj2(l_shell,l,s,j,jm_shell,jm_shell)
               else if (n1 == 1 .and. n2 == 1) then
                  wa = coeflsjj2(l_shell,l,s,j,jm_shell,jp_shell)
               else if (n1 == 0 .and. n2 == 2) then
                  wa = coeflsjj2(l_shell,l,s,j,jp_shell,jp_shell)
               end if
            end if
         else if (l_shell==1 .or. l_shell==2 .or. l_shell==3)  then
            wa = coeflsjjs(l_shell,nn,w,q,l,s,j,n1,qm,jm,qp,jp)
         end if
      else if (l == 0)  then
         if (s == j .and. jm_shell == 1 .and. nn == n1)  then
            wa = one
         else if (s == j .and. jp_shell == 1 .and. nn == n2)  then
            wa = one
         else
            wa = zero
         end if
      else
         stop "coefLSjj: program stop A."
      end if
      wa = factor * wa
      end function coeflsjj
!
!***********************************************************************
!                                                                      *
      recursive function coeflsjj2(l_shell,l,s,j,jm_shell,jp_shell)   result(wa)
!                                                                      *
!     Returns the value of the LS-jj transformation matrix             *
!                                                (l^2 LSJ| j_1 j_2 J). *
!                                                                      *
!     Note that all (generalized) angular momentum quantum numbers     *
!     except for l must be given as twice the original numbers, i.e.   *
!     for the quantum numbers L, S, J, jm_shell, jp_shell.             *
!                                                                      *
!     Calls: nine.                                                     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use numeric_constants,          only: zero, one
      use wigner_9j, only: nine
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: l_shell, l, s, j, jm_shell, jp_shell
      real(real64)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer      :: delta_j
      real(real64) :: rac9
!-----------------------------------------------
      wa = zero
      if (mod(l+s,4) /= 0)  then
         wa = zero
      else if (mod(l+s+j,2) /= 0)  then
         wa = zero
      else if (mod(jm_shell+jp_shell+j,2) /= 0)  then
         wa = zero
      else
         if (jm_shell == jp_shell) then
            if (mod(j,4) /= 0)  then
               wa = zero
            else
               wa = (jm_shell+one) * sqrt((l+one) * (s+one))
            end if
         else
            wa = sqrt(2*(l+one)*(s+one)*(jm_shell+one)*(jp_shell+one))
         end if
         call nine(2*l_shell,2*l_shell,l,1,1,s,jm_shell,jp_shell,j,   &
                   1,delta_j,rac9)
         if (delta_j /= 0)  then
            call nine(2*l_shell,2*l_shell,l,1,1,s,jm_shell,jp_shell,j,&
                      0,delta_j,rac9)
            wa = wa * rac9
         end if
      end if
      end function coeflsjj2
!
!***********************************************************************
!                                                                      *
      recursive function coeflsjjs(lshell,n,w,q,l,s,j,nm,qm,jm,qp,jp)  result(wa)
!                                                                      *
!     Returns the value of the LS-jj transformation matrix for a given *
!     set of quantum numbers.                                          *
!                                                                      *
!     Note that all (generalized) angular momentum quantum numbers     *
!     except for l must be given as twice the original numbers, i.e.   *
!     for the quantum numbers Q, L, S, J, Qm, Jm, Qp, Jp.              *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use ls_jj_subshell_terms, only: ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, ls_jj_number_d5, ls_jj_number_d6, &
         & ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, ls_jj_number_f4, ls_jj_number_f5, &
         & ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, ls_jj_number_p5, ls_jj_number_p6
      use jj2lsj_data_1_c, only: ls_jj_d_10, ls_jj_d_3, ls_jj_d_4, ls_jj_d_5, ls_jj_d_6, ls_jj_d_7, ls_jj_d_8, &
         & ls_jj_d_9, ls_jj_f_3, ls_jj_f_4, ls_jj_f_5, ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, &
         & ls_jj_number_d5, ls_jj_number_d6, ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, &
         & ls_jj_number_f4, ls_jj_number_f5, ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, &
         & ls_jj_number_p5, ls_jj_number_p6, ls_jj_p_3, ls_jj_p_4, ls_jj_p_5, ls_jj_p_6
      use jj2lsj_data_2_c, only: ls_jj_f_6, ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, ls_jj_number_d5, &
         & ls_jj_number_d6, ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, ls_jj_number_f4, &
         & ls_jj_number_f5, ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, ls_jj_number_p5, &
         & ls_jj_number_p6
      use jj2lsj_data_3_c, only: ls_jj_f_7, ls_jj_number_d10, ls_jj_number_d3, ls_jj_number_d4, ls_jj_number_d5, &
         & ls_jj_number_d6, ls_jj_number_d7, ls_jj_number_d8, ls_jj_number_d9, ls_jj_number_f3, ls_jj_number_f4, &
         & ls_jj_number_f5, ls_jj_number_f6, ls_jj_number_f7, ls_jj_number_p3, ls_jj_number_p4, ls_jj_number_p5, &
         & ls_jj_number_p6
      use numeric_constants,          only: zero, one
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: lshell, n, w, q, l, s, j, nm, qm, jm, qp, jp
      real(real64)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      wa = zero
      if (lshell == 0) then
      else if (lshell == 1) then
         select case(n)
         case(3)
!     Use data from the array LS_jj_p_3
            do  i = 1,ls_jj_number_p3
               if(w  ==ls_jj_p_3(i)%w  .and. q  ==ls_jj_p_3(i)%q  .and. &
                  l  ==ls_jj_p_3(i)%l  .and. s  ==ls_jj_p_3(i)%s  .and. &
                  j  ==ls_jj_p_3(i)%j  .and. nm ==ls_jj_p_3(i)%nm .and. &
                  qm ==ls_jj_p_3(i)%qm .and. jm ==ls_jj_p_3(i)%jm .and. &
                  qp ==ls_jj_p_3(i)%qp .and. jp ==ls_jj_p_3(i)%jp) then
                  wa = ls_jj_p_3(i)%factor * &
                       sqrt( one*ls_jj_p_3(i)%nom/ls_jj_p_3(i)%denom )
                  return
               end if
            end do
         case(4)
!     Use data from the array LS_jj_p_4
            do  i = 1,ls_jj_number_p4
               if(w  ==ls_jj_p_4(i)%w  .and. q  ==ls_jj_p_4(i)%q  .and. &
                  l  ==ls_jj_p_4(i)%l  .and. s  ==ls_jj_p_4(i)%s  .and. &
                  j  ==ls_jj_p_4(i)%j  .and. nm ==ls_jj_p_4(i)%nm .and. &
                  qm ==ls_jj_p_4(i)%qm .and. jm ==ls_jj_p_4(i)%jm .and. &
                  qp ==ls_jj_p_4(i)%qp .and. jp ==ls_jj_p_4(i)%jp) then
                  wa = ls_jj_p_4(i)%factor * &
                       sqrt( one*ls_jj_p_4(i)%nom/ls_jj_p_4(i)%denom )
                  return
               end if
            end do
         case(5)
!     Use data from the array LS_jj_p_5
            do  i = 1,ls_jj_number_p5
               if(w  ==ls_jj_p_5(i)%w  .and. q  ==ls_jj_p_5(i)%q  .and. &
                  l  ==ls_jj_p_5(i)%l  .and. s  ==ls_jj_p_5(i)%s  .and. &
                  j  ==ls_jj_p_5(i)%j  .and. nm ==ls_jj_p_5(i)%nm .and. &
                  qm ==ls_jj_p_5(i)%qm .and. jm ==ls_jj_p_5(i)%jm .and. &
                  qp ==ls_jj_p_5(i)%qp .and. jp ==ls_jj_p_5(i)%jp) then
                  wa = ls_jj_p_5(i)%factor * &
                       sqrt( one*ls_jj_p_5(i)%nom/ls_jj_p_5(i)%denom )
                  return
               end if
            end do
         case(6)
!     Use data from the array LS_jj_p_6
            do  i = 1,ls_jj_number_p6
               if(w  ==ls_jj_p_6(i)%w  .and. q  ==ls_jj_p_6(i)%q  .and. &
                  l  ==ls_jj_p_6(i)%l  .and. s  ==ls_jj_p_6(i)%s  .and. &
                  j  ==ls_jj_p_6(i)%j  .and. nm ==ls_jj_p_6(i)%nm .and. &
                  qm ==ls_jj_p_6(i)%qm .and. jm ==ls_jj_p_6(i)%jm .and. &
                  qp ==ls_jj_p_6(i)%qp .and. jp ==ls_jj_p_6(i)%jp) then
                  wa = ls_jj_p_6(i)%factor * &
                       sqrt( one*ls_jj_p_6(i)%nom/ls_jj_p_6(i)%denom )
                  return
               end if
            end do
         case default
            stop "coefLSjjs: program stop A."
         end select
      else if (lshell == 2) then
         select case(n)
         case(3)
!     Use data from the array LS_jj_d_3
            do  i = 1,ls_jj_number_d3
               if(w  ==ls_jj_d_3(i)%w  .and. q  ==ls_jj_d_3(i)%q  .and. &
                  l  ==ls_jj_d_3(i)%l  .and. s  ==ls_jj_d_3(i)%s  .and. &
                  j  ==ls_jj_d_3(i)%j  .and. nm ==ls_jj_d_3(i)%nm .and. &
                  qm ==ls_jj_d_3(i)%qm .and. jm ==ls_jj_d_3(i)%jm .and. &
                  qp ==ls_jj_d_3(i)%qp .and. jp ==ls_jj_d_3(i)%jp) then
                  wa = ls_jj_d_3(i)%factor * &
                       sqrt( one*ls_jj_d_3(i)%nom/ls_jj_d_3(i)%denom )
                  return
               end if
            end do
         case(4)
!     Use data from the array LS_jj_d_4
            do  i = 1,ls_jj_number_d4
               if(w  ==ls_jj_d_4(i)%w  .and. q  ==ls_jj_d_4(i)%q  .and. &
                  l  ==ls_jj_d_4(i)%l  .and. s  ==ls_jj_d_4(i)%s  .and. &
                  j  ==ls_jj_d_4(i)%j  .and. nm ==ls_jj_d_4(i)%nm .and. &
                  qm ==ls_jj_d_4(i)%qm .and. jm ==ls_jj_d_4(i)%jm .and. &
                  qp ==ls_jj_d_4(i)%qp .and. jp ==ls_jj_d_4(i)%jp) then
                  wa = ls_jj_d_4(i)%factor * &
                       sqrt( one*ls_jj_d_4(i)%nom/ls_jj_d_4(i)%denom )
                  return
               end if
            end do
         case(5)
!     Use data from the array LS_jj_d_5
            do  i = 1,ls_jj_number_d5
               if(w  ==ls_jj_d_5(i)%w  .and. q  ==ls_jj_d_5(i)%q  .and. &
                  l  ==ls_jj_d_5(i)%l  .and. s  ==ls_jj_d_5(i)%s  .and. &
                  j  ==ls_jj_d_5(i)%j  .and. nm ==ls_jj_d_5(i)%nm .and. &
                  qm ==ls_jj_d_5(i)%qm .and. jm ==ls_jj_d_5(i)%jm .and. &
                  qp ==ls_jj_d_5(i)%qp .and. jp ==ls_jj_d_5(i)%jp) then
                  wa = ls_jj_d_5(i)%factor * &
                       sqrt( one*ls_jj_d_5(i)%nom/ls_jj_d_5(i)%denom )
                  return
               end if
            end do
         case(6)
!     Use data from the array LS_jj_d_6
            do  i = 1,ls_jj_number_d6
               if(w  ==ls_jj_d_6(i)%w  .and. q  ==ls_jj_d_6(i)%q  .and. &
                  l  ==ls_jj_d_6(i)%l  .and. s  ==ls_jj_d_6(i)%s  .and. &
                  j  ==ls_jj_d_6(i)%j  .and. nm ==ls_jj_d_6(i)%nm .and. &
                  qm ==ls_jj_d_6(i)%qm .and. jm ==ls_jj_d_6(i)%jm .and. &
                  qp ==ls_jj_d_6(i)%qp .and. jp ==ls_jj_d_6(i)%jp) then
                  wa = ls_jj_d_6(i)%factor * &
                       sqrt( one*ls_jj_d_6(i)%nom/ls_jj_d_6(i)%denom )
                  return
               end if
            end do
         case(7)
!     Use data from the array LS_jj_d_7
            do  i = 1,ls_jj_number_d7
               if(w  ==ls_jj_d_7(i)%w  .and. q  ==ls_jj_d_7(i)%q  .and. &
                  l  ==ls_jj_d_7(i)%l  .and. s  ==ls_jj_d_7(i)%s  .and. &
                  j  ==ls_jj_d_7(i)%j  .and. nm ==ls_jj_d_7(i)%nm .and. &
                  qm ==ls_jj_d_7(i)%qm .and. jm ==ls_jj_d_7(i)%jm .and. &
                  qp ==ls_jj_d_7(i)%qp .and. jp ==ls_jj_d_7(i)%jp) then
                  wa = ls_jj_d_7(i)%factor * &
                       sqrt( one*ls_jj_d_7(i)%nom/ls_jj_d_7(i)%denom )
                  return
               end if
            end do
         case(8)
!     Use data from the array LS_jj_d_8
            do  i = 1,ls_jj_number_d8
               if(w  ==ls_jj_d_8(i)%w  .and. q  ==ls_jj_d_8(i)%q  .and. &
                  l  ==ls_jj_d_8(i)%l  .and. s  ==ls_jj_d_8(i)%s  .and. &
                  j  ==ls_jj_d_8(i)%j  .and. nm ==ls_jj_d_8(i)%nm .and. &
                  qm ==ls_jj_d_8(i)%qm .and. jm ==ls_jj_d_8(i)%jm .and. &
                  qp ==ls_jj_d_8(i)%qp .and. jp ==ls_jj_d_8(i)%jp) then
                  wa = ls_jj_d_8(i)%factor * &
                       sqrt( one*ls_jj_d_8(i)%nom/ls_jj_d_8(i)%denom )
                  return
               end if
            end do
         case(9)
!     Use data from the array LS_jj_d_9
            do  i = 1,ls_jj_number_d9
               if(w  ==ls_jj_d_9(i)%w  .and. q  ==ls_jj_d_9(i)%q  .and. &
                  l  ==ls_jj_d_9(i)%l  .and. s  ==ls_jj_d_9(i)%s  .and. &
                  j  ==ls_jj_d_9(i)%j  .and. nm ==ls_jj_d_9(i)%nm .and. &
                  qm ==ls_jj_d_9(i)%qm .and. jm ==ls_jj_d_9(i)%jm .and. &
                  qp ==ls_jj_d_9(i)%qp .and. jp ==ls_jj_d_9(i)%jp) then
                  wa = ls_jj_d_9(i)%factor * &
                       sqrt( one*ls_jj_d_9(i)%nom/ls_jj_d_9(i)%denom )
                  return
               end if
            end do
         case(10)
!     Use data from the array LS_jj_d_10
            do  i = 1,ls_jj_number_d10
               if(w ==ls_jj_d_10(i)%w  .and. q ==ls_jj_d_10(i)%q  .and. &
                  l ==ls_jj_d_10(i)%l  .and. s ==ls_jj_d_10(i)%s  .and. &
                  j ==ls_jj_d_10(i)%j  .and. nm==ls_jj_d_10(i)%nm .and. &
                  qm==ls_jj_d_10(i)%qm .and. jm==ls_jj_d_10(i)%jm .and. &
                  qp==ls_jj_d_10(i)%qp .and. jp==ls_jj_d_10(i)%jp) then
                  wa= ls_jj_d_10(i)%factor * &
                       sqrt( one*ls_jj_d_10(i)%nom/ls_jj_d_10(i)%denom )
                  return
               end if
            end do
         case default
            stop "coefLSjjs: program stop B."
         end select
      else if (lshell == 3) then
         select case(n)
         case(3)
!     Use data from the array LS_jj_f_3
            do  i = 1,ls_jj_number_f3
               if(w  ==ls_jj_f_3(i)%w  .and. q  ==ls_jj_f_3(i)%q  .and. &
                  l  ==ls_jj_f_3(i)%l  .and. s  ==ls_jj_f_3(i)%s  .and. &
                  j  ==ls_jj_f_3(i)%j  .and. nm ==ls_jj_f_3(i)%nm .and. &
                  qm ==ls_jj_f_3(i)%qm .and. jm ==ls_jj_f_3(i)%jm .and. &
                  qp ==ls_jj_f_3(i)%qp .and. jp ==ls_jj_f_3(i)%jp) then
                  wa = ls_jj_f_3(i)%factor * &
                       sqrt( one*ls_jj_f_3(i)%nom/ls_jj_f_3(i)%denom )
                  return
               end if
            end do
         case(4)
!     Use data from the array LS_jj_f_4
            do  i = 1,ls_jj_number_f4
               if(w  ==ls_jj_f_4(i)%w  .and. q  ==ls_jj_f_4(i)%q  .and. &
                  l  ==ls_jj_f_4(i)%l  .and. s  ==ls_jj_f_4(i)%s  .and. &
                  j  ==ls_jj_f_4(i)%j  .and. nm ==ls_jj_f_4(i)%nm .and. &
                  qm ==ls_jj_f_4(i)%qm .and. jm ==ls_jj_f_4(i)%jm .and. &
                  qp ==ls_jj_f_4(i)%qp .and. jp ==ls_jj_f_4(i)%jp) then
                  wa = ls_jj_f_4(i)%factor * &
                       sqrt( one*ls_jj_f_4(i)%nom/ls_jj_f_4(i)%denom )
                  return
               end if
            end do
         case(5)
!     Use data from the array LS_jj_f_5
            do  i = 1,ls_jj_number_f5
               if(w  ==ls_jj_f_5(i)%w  .and. q  ==ls_jj_f_5(i)%q  .and. &
                  l  ==ls_jj_f_5(i)%l  .and. s  ==ls_jj_f_5(i)%s  .and. &
                  j  ==ls_jj_f_5(i)%j  .and. nm ==ls_jj_f_5(i)%nm .and. &
                  qm ==ls_jj_f_5(i)%qm .and. jm ==ls_jj_f_5(i)%jm .and. &
                  qp ==ls_jj_f_5(i)%qp .and. jp ==ls_jj_f_5(i)%jp) then
                  wa = ls_jj_f_5(i)%factor * &
                       sqrt( one*ls_jj_f_5(i)%nom/ls_jj_f_5(i)%denom )
                  return
               end if
            end do
         case(6)
!     Use data from the array LS_jj_f_6
            do  i = 1,ls_jj_number_f6
               if(w  ==ls_jj_f_6(i)%w  .and. q  ==ls_jj_f_6(i)%q  .and. &
                  l  ==ls_jj_f_6(i)%l  .and. s  ==ls_jj_f_6(i)%s  .and. &
                  j  ==ls_jj_f_6(i)%j  .and. nm ==ls_jj_f_6(i)%nm .and. &
                  qm ==ls_jj_f_6(i)%qm .and. jm ==ls_jj_f_6(i)%jm .and. &
                  qp ==ls_jj_f_6(i)%qp .and. jp ==ls_jj_f_6(i)%jp) then
                  wa = ls_jj_f_6(i)%factor * &
                       sqrt( one*ls_jj_f_6(i)%nom/ls_jj_f_6(i)%denom )
                  return
               end if
            end do
         case(7)
!     Use data from the array LS_jj_f_7
            do  i = 1, ls_jj_number_f7, 1
               if (w == ls_jj_f_7(i)%w) then
                  if(q == ls_jj_f_7(i)%q) then
                     if(l == ls_jj_f_7(i)%l) then
                        if(s == ls_jj_f_7(i)%s) then
                           if(j == ls_jj_f_7(i)%j) then
                              if(nm == ls_jj_f_7(i)%nm) then
                                 if(qm == ls_jj_f_7(i)%qm .and. &
                                    jm == ls_jj_f_7(i)%jm .and. &
                                    qp == ls_jj_f_7(i)%qp .and. &
                                    jp == ls_jj_f_7(i)%jp) then
                                       wa = one * ls_jj_f_7(i)%factor  * &
                                            dsqrt(one*ls_jj_f_7(i)%nom)/ &
                                            dsqrt(one*ls_jj_f_7(i)%denom )
                                       return
                                 end if
                              end if
                           end if
                        end if
                     end if
                  end if
               end if
            end do
         case default
            stop "coefLSjjs: program stop C."
         end select
      else
         stop "coefLSjjs: program stop D."
      end if
      end function coeflsjjs
!
!***********************************************************************
!                                                                      *
      subroutine dallocasfls(asf_set_ls)
!                                                                      *
!     Dellocates the storage of asf_set_LS.                            *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use selected_level_indices,             only: nvec
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: itjpo
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      type(asf_basis_ls), intent(inout) :: asf_set_ls
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: i
!-----------------------------------------------
      deallocate(asf_set_ls%csf_set_ls%shell)
      deallocate(asf_set_ls%csf_set_ls%parent)
      do i = 1, asf_set_ls%csf_set_ls%nocsf, 1
         deallocate(asf_set_ls%csf_set_ls%csf(i)%occupation)
         deallocate(asf_set_ls%csf_set_ls%csf(i)%seniority)
         deallocate(asf_set_ls%csf_set_ls%csf(i)%w)
         deallocate(asf_set_ls%csf_set_ls%csf(i)%shelll)
         deallocate(asf_set_ls%csf_set_ls%csf(i)%shells)
         deallocate(asf_set_ls%csf_set_ls%csf(i)%shelllx)
         deallocate(asf_set_ls%csf_set_ls%csf(i)%shellsx)
      end do
      deallocate(asf_set_ls%csf_set_ls%csf)
      do i = 1, nvec
         deallocate(asf_set_ls%asf(i)%eigenvector)
      end do
      deallocate(asf_set_ls%asf)
      end subroutine dallocasfls
!
!***********************************************************************
!                                                                      *
      subroutine getchls(csf_index,shell_number,string_ls, string_xls)
!                                                                      *
!     The spectroscopic notation of a shell in LS coupling is returned.*
!                                                                      *
!     Calls: convrt_double.                                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      use half_integer_format, only: convrt_double
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: csf_index,shell_number
      character(len=4), intent(out) :: string_ls, string_xls
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: text_length
      character(len=1) :: string_s, string_v
      character(len=4) :: string_cnum
!-----------------------------------------------
      call                                                            &
      convrt_double(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%shells(shell_number)),&
      string_cnum,text_length)
      string_s = string_cnum(1:text_length)
      if (asf_set_ls%csf_set_ls%shell(shell_number)%l < 3)    then
         call                                                         &
         convrt_double(2*asf_set_ls%csf_set_ls%csf(csf_index)%seniority(shell_number),&
         string_cnum,text_length)
         string_v = string_cnum(1:text_length)
      else
         call                                                         &
         convrt_double(2*asf_set_ls%csf_set_ls%csf(csf_index)%w(shell_number),&
         string_cnum,text_length)
         string_v = string_cnum(1:text_length)
      end if
      string_ls = string_s //                                         &
      l_string(asf_set_ls%csf_set_ls%csf(csf_index)%shelll(shell_number)/2)&
      //string_v
      !
      call                                                            &
      convrt_double(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%shellsx(shell_number)),&
      string_cnum,text_length)
      string_s =  string_cnum(1:text_length)
      string_xls = string_s //                                        &
         l_string(asf_set_ls%csf_set_ls%csf(csf_index)%shelllx(shell_number)/2)
      end subroutine getchls
!
!***********************************************************************
!                                                                      *
      subroutine getchls_lbl(csf_index,shell_number,string_ls,string_xls)
!                                                                      *
!     The spectroscopic notation of a shell in LS coupling is returned.*
!                                                                      *
!     Calls: convrt_double.                                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     Modified by G. Gaigalas                                   2024   *
!                                                                      *
!***********************************************************************
      use half_integer_format, only: convrt_double
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)           :: csf_index,shell_number
      character(len=4), intent(out) :: string_ls, string_xls
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer          :: text_length, w, tl, ts, seniority
      character(len=1) :: string_s, string_v
      character(len=4) :: string_cnum
!-----------------------------------------------
      call                                                            &
      convrt_double(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%shells(shell_number)),&
      string_cnum,text_length)
      string_s = string_cnum(1:text_length)
!GG-2024
      if (asf_set_ls%csf_set_ls%shell(shell_number)%l == 2)    then
         if (ls_convention == 0) then
            call                                                       &
                convrt_double(2*asf_set_ls%csf_set_ls%csf(csf_index)% &
                seniority(shell_number),string_cnum,text_length)
            string_v = string_cnum(1:text_length)
         else
           ts=asf_set_ls%csf_set_ls%csf(csf_index)%shells(shell_number)
           tl=asf_set_ls%csf_set_ls%csf(csf_index)%shelll(shell_number)
           seniority=asf_set_ls%csf_set_ls%csf(csf_index)%seniority(shell_number)
           if(ts+1 == 1) then
              if(tl == 0) then
! 1S
                 if(seniority == 0) then
                    w = 1
                 else if(seniority == 4) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 1S."
                 end if
              else if(tl == 4) then
! 1D
                 if(seniority == 2) then
                    w = 1
                 else if(seniority == 4) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 1D."
                 end if
              else if(tl == 8) then
! 1G
                 if(seniority == 2) then
                    w = 1
                 else if(seniority == 4) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 1G."
                 end if
              else
                 w = 0
              end if
           else if(ts+1 == 2) then
              if(tl == 4) then
! 2D
                 if(seniority == 1) then
                    w = 1
                 else if(seniority == 3) then
                    w = 2
                 else if(seniority == 5) then
                    w = 3
                 else
                    stop "getchLS_lbl(): program stop 2D."
                 end if
              else if(tl == 6) then
! 2F
                 if(seniority == 3) then
                    w = 1
                 else if(seniority == 5) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 2F."
                 end if
              else if(tl == 8) then
! 2G
                 if(seniority == 3) then
                    w = 1
                 else if(seniority == 5) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 2G."
                 end if
              else
                 w = 0
              end if
           else if(ts+1 == 3) then
              if(tl == 2) then
! 3P
                 if(seniority == 2) then
                    w = 1
                 else if(seniority == 4) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 3P."
                 end if
              else if(tl == 6) then
! 3F
                 if(seniority == 2) then
                    w = 1
                 else if(seniority == 4) then
                    w = 2
                 else
                    stop "getchLS_lbl(): program stop 3F."
                 end if
              else
                 w = 0
              end if
           else
              w = 0
           end if
           call convrt_double(2*w,string_cnum,text_length)
           string_v = string_cnum(1:text_length)
         end if
      else if (asf_set_ls%csf_set_ls%shell(shell_number)%l < 3)    then
         if (ls_convention == 0) then
            call                                                       &
                convrt_double(2*asf_set_ls%csf_set_ls%csf(csf_index)% &
                seniority(shell_number),string_cnum,text_length)
            string_v = string_cnum(1:text_length)
         else
            w = asf_set_ls%csf_set_ls%csf(csf_index)%w(shell_number)
            call convrt_double(2*w,string_cnum,text_length)
            string_v = string_cnum(1:text_length)
         end if
      else
         w = asf_set_ls%csf_set_ls%csf(csf_index)%w(shell_number)
         call convrt_double(2*w, string_cnum,text_length)
         string_v = string_cnum(1:text_length)
      end if
      if(asf_set_ls%csf_set_ls%shell(shell_number)%l < 2) then
         string_ls = string_s //                                      &
         l_string(asf_set_ls%csf_set_ls%csf(csf_index)%              &
                                              shelll(shell_number)/2)
      else if(asf_set_ls%csf_set_ls%shell(shell_number)%l == 3        &
                                                   .and. w == 0) then
         string_ls = string_s //                                      &
         l_string(asf_set_ls%csf_set_ls%csf(csf_index)%              &
                                              shelll(shell_number)/2)
      else if(ls_convention == 1 .and. w == 0) then
         string_ls = string_s //                                      &
         l_string(asf_set_ls%csf_set_ls%csf(csf_index)%              &
                                              shelll(shell_number)/2)
      else
         if(w == 10) then
            string_ls = string_s //                                   &
            l_string(asf_set_ls%csf_set_ls%csf(csf_index)%shelll(    &
            shell_number)/2)//"A"
         else
            string_ls = string_s //                                   &
            l_string(asf_set_ls%csf_set_ls%csf(csf_index)%shelll(    &
            shell_number)/2)//string_v
         end if
      end if
!GG-2024 end
      !
      call                                                            &
      convrt_double(2*(1+asf_set_ls%csf_set_ls%csf(csf_index)%shellsx(shell_number)),&
      string_cnum,text_length)
      string_s =  string_cnum(1:text_length)
      string_xls = string_s //                                        &
         l_string(asf_set_ls%csf_set_ls%csf(csf_index)%shelllx(shell_number)/2)
      end subroutine getchls_lbl
!
!***********************************************************************
!                                                                      *
      subroutine getxj
!                                                                      *
!                                                                      *
!     Calls: jcup, jqs, ichop.                                         *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use csf_pair_angular_state,             only: core_orbital_count => ncore
      use orbital_csf_basis,           only: csf_count, orbital_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: jqs, jcup
      use subshell_occupation_class, only: ichop
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: jcnt, jcntop, csf_index, orbital_index
!-----------------------------------------------
      jcoup = 0
      do csf_index = 1, csf_count
         jcnt = 1
         jcntop = 0
         do orbital_index = core_orbital_count+1, orbital_count
            if(orbital_index == 1) then
               if(ichop(orbital_index,csf_index) /= 0) then
                  if(orbital_count == 1) then
                      jcoup(orbital_index,csf_index) = 0
                  else
                      if(ichop(orbital_index+1,csf_index) /= 0) then
                         jcoup(orbital_index,csf_index)   = 0
                         jcoup(orbital_index+1,csf_index) = 0
                      else
                         jcoup(orbital_index,csf_index)   = 0
                         jcoup(orbital_index+1,csf_index) = jqs(3,orbital_index+1,csf_index) - 1
                         jcntop = 1
                     end if
                  end if
               else
                  jcntop = 1
                  if(orbital_count > 1) then
                     if (ichop(orbital_index+1,csf_index) == 0) then
!GG 2015_06_21 Gediminas Gaigalas
!GG                        Jcoup(JNW,JNCF)   = JCUP(JCNT,JNCF) - 1
!CYC-2023
! Modifications for configurations with empty 1s subshell
!                        Jcoup(JNW,JNCF)   = JQS(3,JNW,JNCF) - 1
                        jcoup(orbital_index,csf_index)   = iabs(jcup(jcnt,csf_index) - 2)
!CYC-2023 end
                        jcoup(orbital_index+1,csf_index) = jcup(jcnt,csf_index) - 1
                        jcnt = jcnt + 1
                     else
                        jcoup(orbital_index,csf_index)   = jqs(3,orbital_index,csf_index) - 1
                        jcoup(orbital_index+1,csf_index) = jqs(3,orbital_index,csf_index) - 1
                     end if
                  else
                     jcoup(orbital_index,csf_index) = jcup(jcnt,csf_index) - 1
                  end if
               end if
            else if(orbital_index == 2 .and. core_orbital_count+1 == 2) then
               if(ichop(orbital_index,csf_index) /= 0) then
                  jcoup(orbital_index,csf_index) = 0
               else
                  jcntop = 1
                  jcoup(orbital_index,csf_index)   = jqs(3,orbital_index,csf_index) - 1
               end if
            else if(orbital_index > 2) then
               if (ichop(orbital_index,csf_index) /= 0) then
                  if(orbital_index == core_orbital_count+1) then
                     jcoup(orbital_index,csf_index) = jqs(3,orbital_index,csf_index) - 1
                  else
                     jcoup(orbital_index,csf_index) = jcoup(orbital_index-1,csf_index)
                  end if
               else
                  if (jcntop /= 0) then
                     jcoup(orbital_index,csf_index) = jcup(jcnt,csf_index)  - 1
                     jcnt = jcnt + 1
                  else
                     jcoup(orbital_index,csf_index) = jqs(3,orbital_index,csf_index) - 1
                  end if
                  jcntop = jcntop + 1
               end if
            end if
         end do
      end do
      return
      end subroutine  getxj
!
!***********************************************************************
!                                                                      *
      subroutine gettermls(l_shell,n,ls,number)
!                                                                      *
!     This procedure returns all allowed subshell terms (l,w,Q,L,S)    *
!     for given l^N, N = 0, 1, 2 or 3.                                 *
!                                                                      *
!     Calls:                                                           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use ls_jj_subshell_terms, only: subshell_term_ls, term_ls_d, term_ls_f, term_ls_g1, term_ls_g2, term_ls_h1, term_ls_h2, &
         & term_ls_i1, term_ls_i2, term_ls_k1, term_ls_k2, term_ls_l1, term_ls_l2, term_ls_m1, term_ls_m2, &
         & term_ls_p, term_ls_s
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)                                 :: l_shell, n
      type(subshell_term_ls), intent(out) :: ls(120)
      integer, intent(out)                                :: number
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: m_q, i, j
!-----------------------------------------------
      m_q = n - 2* l_shell - 1;  j = 0
      select case (l_shell)
      case (0)
         do i = 1,2
            if (mod(m_q + term_ls_s(i)%subshell_two_quasispin,2) == 0) then
               if (abs(m_q) <= term_ls_s(i)%subshell_two_quasispin) then
                  j = j + 1
                  ls(j)%subshell_l  = term_ls_s(i)%subshell_l
                  ls(j)%w        = term_ls_s(i)%w
                  ls(j)%subshell_two_quasispin        = term_ls_s(i)%subshell_two_quasispin
                  ls(j)%subshell_two_l       = term_ls_s(i)%subshell_two_l
                  ls(j)%subshell_two_s        = term_ls_s(i)%subshell_two_s
               end if
            end if
         end do
      case (1)
         do i = 1,6
            if (mod(m_q + term_ls_p(i)%subshell_two_quasispin,2) == 0) then
               if (abs(m_q) <= term_ls_p(i)%subshell_two_quasispin) then
                  j = j + 1
                  ls(j)%subshell_l  = term_ls_p(i)%subshell_l
                  ls(j)%w        = term_ls_p(i)%w
                  ls(j)%subshell_two_quasispin        = term_ls_p(i)%subshell_two_quasispin
                  ls(j)%subshell_two_l       = term_ls_p(i)%subshell_two_l
                  ls(j)%subshell_two_s        = term_ls_p(i)%subshell_two_s
               end if
            end if
         end do
      case (2)
         do i = 1,32
            if (mod(m_q + term_ls_d(i)%subshell_two_quasispin,2) == 0) then
               if (abs(m_q) <= term_ls_d(i)%subshell_two_quasispin) then
                  j = j + 1
                  ls(j)%subshell_l  = term_ls_d(i)%subshell_l
                  ls(j)%w        = term_ls_d(i)%w
                  ls(j)%subshell_two_quasispin        = term_ls_d(i)%subshell_two_quasispin
                  ls(j)%subshell_two_l       = term_ls_d(i)%subshell_two_l
                  ls(j)%subshell_two_s        = term_ls_d(i)%subshell_two_s
               end if
            end if
         end do
      case (3)
         do i = 1,238
            if (mod(m_q + term_ls_f(i)%subshell_two_quasispin,2) == 0) then
               if (abs(m_q) <= term_ls_f(i)%subshell_two_quasispin) then
                  j = j + 1
                  ls(j)%subshell_l  = term_ls_f(i)%subshell_l
                  ls(j)%w        = term_ls_f(i)%w
                  ls(j)%subshell_two_quasispin        = term_ls_f(i)%subshell_two_quasispin
                  ls(j)%subshell_two_l       = term_ls_f(i)%subshell_two_l
                  ls(j)%subshell_two_s        = term_ls_f(i)%subshell_two_s
               end if
            end if
         end do
      case (4)
         select case (n)
         case (1)
            i = 1;       j = 1
            ls(j)%subshell_l  = term_ls_g1(i)%subshell_l
            ls(j)%w        = term_ls_g1(i)%w
            ls(j)%subshell_two_quasispin        = term_ls_g1(i)%subshell_two_quasispin
            ls(j)%subshell_two_l       = term_ls_g1(i)%subshell_two_l
            ls(j)%subshell_two_s        = term_ls_g1(i)%subshell_two_s
         case (2)
            do i = 1,9
               j = j + 1
               ls(j)%subshell_l  = term_ls_g2(i)%subshell_l
               ls(j)%w        = term_ls_g2(i)%w
               ls(j)%subshell_two_quasispin        = term_ls_g2(i)%subshell_two_quasispin
               ls(j)%subshell_two_l       = term_ls_g2(i)%subshell_two_l
               ls(j)%subshell_two_s        = term_ls_g2(i)%subshell_two_s
            end do
         case default
            stop  "gettermLS(): program stop A."
         end select
      case (5)
         select case (n)
         case (1)
            i = 1;       j = 1
            ls(j)%subshell_l  = term_ls_h1(i)%subshell_l
            ls(j)%w        = term_ls_h1(i)%w
            ls(j)%subshell_two_quasispin        = term_ls_h1(i)%subshell_two_quasispin
            ls(j)%subshell_two_l       = term_ls_h1(i)%subshell_two_l
            ls(j)%subshell_two_s        = term_ls_h1(i)%subshell_two_s
         case (2)
            do i = 1,11
               j = j + 1
               ls(j)%subshell_l  = term_ls_h2(i)%subshell_l
               ls(j)%w        = term_ls_h2(i)%w
               ls(j)%subshell_two_quasispin        = term_ls_h2(i)%subshell_two_quasispin
               ls(j)%subshell_two_l       = term_ls_h2(i)%subshell_two_l
               ls(j)%subshell_two_s        = term_ls_h2(i)%subshell_two_s
            end do
         case default
            stop  "gettermLS(): program stop B."
         end select
      case (6)
         select case (n)
         case (1)
            i = 1;       j = 1
            ls(j)%subshell_l  = term_ls_i1(i)%subshell_l
            ls(j)%w        = term_ls_i1(i)%w
            ls(j)%subshell_two_quasispin        = term_ls_i1(i)%subshell_two_quasispin
            ls(j)%subshell_two_l       = term_ls_i1(i)%subshell_two_l
            ls(j)%subshell_two_s        = term_ls_i1(i)%subshell_two_s
         case (2)
            do i = 1,13
               j = j + 1
               ls(j)%subshell_l  = term_ls_i2(i)%subshell_l
               ls(j)%w        = term_ls_i2(i)%w
               ls(j)%subshell_two_quasispin        = term_ls_i2(i)%subshell_two_quasispin
               ls(j)%subshell_two_l       = term_ls_i2(i)%subshell_two_l
               ls(j)%subshell_two_s        = term_ls_i2(i)%subshell_two_s
            end do
         case default
            stop  "gettermLS(): program stop C."
         end select
      case (7)
         select case (n)
         case (1)
            i = 1;       j = 1
            ls(j)%subshell_l  = term_ls_k1(i)%subshell_l
            ls(j)%w        = term_ls_k1(i)%w
            ls(j)%subshell_two_quasispin        = term_ls_k1(i)%subshell_two_quasispin
            ls(j)%subshell_two_l       = term_ls_k1(i)%subshell_two_l
            ls(j)%subshell_two_s        = term_ls_k1(i)%subshell_two_s
         case (2)
            do i = 1,15
               j = j + 1
               ls(j)%subshell_l  = term_ls_k2(i)%subshell_l
               ls(j)%w        = term_ls_k2(i)%w
               ls(j)%subshell_two_quasispin        = term_ls_k2(i)%subshell_two_quasispin
               ls(j)%subshell_two_l       = term_ls_k2(i)%subshell_two_l
               ls(j)%subshell_two_s        = term_ls_k2(i)%subshell_two_s
            end do
         case default
            stop  "gettermLS(): program stop D."
         end select
      case (8)
         select case (n)
         case (1)
            i = 1;       j = 1
            ls(j)%subshell_l  = term_ls_l1(i)%subshell_l
            ls(j)%w        = term_ls_l1(i)%w
            ls(j)%subshell_two_quasispin        = term_ls_l1(i)%subshell_two_quasispin
            ls(j)%subshell_two_l       = term_ls_l1(i)%subshell_two_l
            ls(j)%subshell_two_s        = term_ls_l1(i)%subshell_two_s
         case (2)
            do i = 1,17
               j = j + 1
               ls(j)%subshell_l  = term_ls_l2(i)%subshell_l
               ls(j)%w        = term_ls_l2(i)%w
               ls(j)%subshell_two_quasispin        = term_ls_l2(i)%subshell_two_quasispin
               ls(j)%subshell_two_l       = term_ls_l2(i)%subshell_two_l
               ls(j)%subshell_two_s        = term_ls_l2(i)%subshell_two_s
            end do
         case default
            stop  "gettermLS(): program stop E."
         end select
      case (9)
         select case (n)
         case (1)
            i = 1;       j = 1
            ls(j)%subshell_l  = term_ls_m1(i)%subshell_l
            ls(j)%w        = term_ls_m1(i)%w
            ls(j)%subshell_two_quasispin        = term_ls_m1(i)%subshell_two_quasispin
            ls(j)%subshell_two_l       = term_ls_m1(i)%subshell_two_l
            ls(j)%subshell_two_s        = term_ls_m1(i)%subshell_two_s
         case (2)
            do i = 1,19
               j = j + 1
               ls(j)%subshell_l  = term_ls_m2(i)%subshell_l
               ls(j)%w        = term_ls_m2(i)%w
               ls(j)%subshell_two_quasispin        = term_ls_m2(i)%subshell_two_quasispin
               ls(j)%subshell_two_l       = term_ls_m2(i)%subshell_two_l
               ls(j)%subshell_two_s        = term_ls_m2(i)%subshell_two_s
            end do
         case default
            stop  "gettermLS(): program stop F."
         end select
      case default
         stop  "gettermLS(): program stop G."
      end select
      number = j
      end subroutine gettermls
!
!***********************************************************************
!                                                                      *
      subroutine inscreen(thresh,levels,number_of_levels,ioutc,ioutj,  &
                                                          unique,ioutt)
      use, intrinsic :: iso_fortran_env, only: real32
!                                                                      *
!     The input from the screen.                                       *
!                                                                      *
!     Calls: sercsla, getxj, getmixblock, inscreenlev, openfl.         *
!     inscreenlev.                                                     *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!     Modified by G. Gaigalas                                   2022   *
!     Modified by G. Gaigalas                                   2024   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis,           only: csf_count, orbital_count
      use selected_level_indices,          only: nvec
      use standard_io_units,        only: stdin_unit, stderr_unit
      use numeric_constants,          only: eps, zero
      use transition_basis_blocks,           only: nevinblk, ncfinblk, nblock, two_j
      use csf_pair_angular_state,             only: core_orbital_count => ncore
      use calculation_state_files, only: read_state_name
      use atomic_calculation_state,           only: nuclear_charge, electron_count
      use getmixblock_i, only: getmixblock
      use yes_no_input, only: getyn
      use file_open, only: openfl
      use csf_basis_loader, only: setcsla
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(out)      :: ioutc,ioutj,unique,ioutt
      real(real64), intent(out) :: thresh
      integer, intent(out) :: number_of_levels(blocks_number)
      integer, intent(out) :: levels(blocks_number, vectors_number)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: phase_start
      integer            :: nelecdum, ncftotdum, nwdum, nblockdum
      integer            :: jbdum, ncfinblkdum, nevinblkdum, iatjpdum
      integer            :: i, ii, isum, k, nci, io_status, jb
      integer            :: iblock, number_of_levels_tmp
      logical            :: yes, fail
      character(len=6)   :: mixing_file_signature
      character(len=24)  :: name
      character(len=256) :: record, util_csl_file
      integer    :: posi(blocks_number)
      integer :: levels_tmp(1:vectors_number)
!-----------------------------------------------
      nblock = 0
      call read_state_name(name)
      k = len_trim(name) + 1
!
!   Open, check, load data from, and close, the  .csl  file
      phase_start=profile_time()
      call setcsla(name,core_orbital_count)
      call profile_emit("csf_read_parse",phase_start,items=csf_count)
      phase_start=profile_time()
      allocate(jcoup(1:orbital_count,1:csf_count))
      call getxj
      call profile_emit("jcoup",phase_start,items=csf_count)
      write (stderr_unit,*)
      write (stderr_unit,*) "Mixing coefficients from a CI calc.?"
      yes = getyn ()
      if (yes) then
         nci = 0
      else
         nci = 1
      end if
!GG-2017
      write (stderr_unit,*) "Do you need a unique labeling? (y/n)"
      yes = getyn ()
      if (yes) then
         unique = 1
      else
         unique = 0
      end if
!GG-2017 end
!
!   Get the eigenvectors
      phase_start=profile_time()
      call getmixblock (name, nci)
      call profile_emit("mix_read",phase_start,items=nvec)
      write (stderr_unit,*) "Default settings?  (y/n) "
      yes = getyn ()
      if (yes) then
         epsnew = 0.005_real64
         isum = 0
         imincompoff = 1
         mincomp = 1
         thresh = 0.001_real64
         ioutc = 0
         ioutj = 0
         ioutt = 0
         do i = 1, nblock
            number_of_levels(i) = nevinblk(i)
            if(nevinblk(i) /= 0) then
               do ii = 1, nevinblk(i)
                  isum = isum + 1
                  levels(i,ii) = isum
               end do
           end if
         end do
!GG-2024
         ls_convention = 1
!GG-2024 end
      else
         write (stderr_unit,*) "All levels (Y/N)"
         yes = getyn ()
         write (stderr_unit,*)
         if (yes) then
            isum = 0
            do i = 1, nblock
                number_of_levels(i) = nevinblk(i)
                if(nevinblk(i) /= 0) then
                   do ii = 1, nevinblk(i)
                      isum = isum + 1
                      levels(i,ii) = isum
                   end do
               end if
            end do
         else
            number_of_levels = 0
            posi(1) = 0
            do i = 1, nblock-1
               posi(i+1) = nevinblk(i) + posi(i)
            end do
            write (stderr_unit,*) "Maximum number of ASFs is:",nvec
            write (stderr_unit,*) "Enter the level numbers of the ASF which are to be transformed,"
    2       write (stderr_unit,*) "Enter the block number"
            read (*, "(I2)") iblock
            write (stderr_unit,*) "The block number is:",iblock
            write (stderr_unit,*) " e.g. 1 3 4  7 - 20  48  69 - 85 :"
            read (*, "(a)") record
            call inscreenlev(record,levels_tmp,number_of_levels_tmp,fail)
            if (fail) then
               write (stderr_unit,*) "Unable to interprete the serial level numbers; redo ..."
               goto 2
            end if
            number_of_levels(iblock) = number_of_levels_tmp
            if (nvec < number_of_levels(iblock)) then
               write (stderr_unit,*) "There are to much ASF:", number_of_levels(iblock)
               goto 2
            end if
            do i = 1,number_of_levels_tmp
               levels(iblock,i) = posi(iblock) + levels_tmp(i)
            end do
            write (stderr_unit,*) ""
            write (stderr_unit,*) "Do you need to include more levels?  (y/n)"
            yes = getyn ()
            if (yes) goto 2
         end if
    3    write (stderr_unit,*) "Maximum % of omitted composition"
!GG-2024
         write (stderr_unit,*) "(default setting: 1.000)"
!GG-2024 end
         read *, mincomp
         if(mincomp == zero) then
            imincompoff = 0
            epsnew = eps*eps
            ioutc = 1
            ioutj = 1
         else if( mincomp > zero) then
            imincompoff = 1
            write (stderr_unit,*) "What is the value below which an eigenvector component"
            write (stderr_unit,*) "is to be neglected in the determination of the LSJ expansion:"
            write (stderr_unit,"(A,F8.5)") " should be smaller than:",mincomp*0.01_real64
!GG-2024
         write (stderr_unit,*) "(default setting: 0.005)"
!GG-2024 end
            read *, epsnew
         else
            write (stderr_unit,*) " The maximum of omitted composition can be 100%"
            goto 3
         end if
         write (stderr_unit,*) "What is the value below which an eigenvector composition"
         write (stderr_unit,*) "is to be neglected for printing?"
!GG-2024
         write (stderr_unit,*) "(default setting: 0.001)"
!GG-2024 end
         read *, thresh
         if( mincomp > zero) then
            write (stderr_unit,*) "Do you need the output file *.lsj.c?  (y/n)"
            yes = getyn ()
            if (yes) then
               ioutc = 1
            else
               ioutc = 0
            end if
            write (stderr_unit,*) "Do you need the output file *.lsj.j?  (y/n)"
            yes = getyn ()
            if (yes) then
               ioutj = 1
            else
               ioutj = 0
            end if
         end if
         if(mincomp == zero) then
            write (stderr_unit,*)                                            &
            "Do you need the transformation output file *.lsj.T?  (y/n)"
            yes = getyn ()
            if (yes) then
               ioutt = 1
            else
               write (stderr_unit,*)                                         &
               "Will you use the transformation file *.lsj.T?  (y/n)"
               yes = getyn ()
               if (yes) then
                  ioutt = 2
               else
                  ioutt = 0
               end if
            end if
         else
            ioutt = 0
         end if
!GG-2024
         write (stderr_unit,*)                                               &
         "Which convention should be used for d shell identification?"
         write (stderr_unit,*)                                               &
           "Seniority notation / Nielson-Koster indexes (0/1)"
         read *, ls_convention
!GG-2024 end
      end if
!
      write (stderr_unit,*)
      if (imincompoff == 1)             &
         write (*,"(A,F8.3)") " Maximum % of omitted composition is ",mincomp
      write (*,"(A,ES8.1,A)")   &
        " Below ",epsnew," the eigenvector component is to be neglected for calculating"
      write (*,"(A,ES8.1,A)")   &
        " Below ",thresh," the eigenvector composition is to be neglected for printing"
!GG-2024
      write (stderr_unit,*)
      if (ls_convention == 0) then
        write (*,"(A)")                                                &
            " Seniority notaton is used for d shell identification in"
        write (*,"(5A)")"   ",name(1:k-1)//".lsj"//".lbl", ", ",       &
                              name(1:k-1)//".uni.lsj"//".lbl",","
        write (*,"(5A)")"   ",name(1:k-1)//".uni.lsj"//".sum",", and ",&
                              name(1:k-1)//".lsj"//".j","  files"
      else
        write (*,"(A)")                                                &
         " Nielson-Koster indexes is used for d shell identification in"
        write (*,"(5A)")"   ",name(1:k-1)//".lsj"//".lbl", ", ",       &
                              name(1:k-1)//".uni.lsj"//".lbl",","
        write (*,"(5A)")"   ",name(1:k-1)//".uni.lsj"//".sum",", and ",&
                              name(1:k-1)//".lsj"//".j","  files"
      end if
      write (*,"(A)")                                                  &
      " while seniority notation for this shell is used in other files"
!GG-2024 end
      print *, " "
!
!     Opening the files  *.lsj.c
!
      if(ioutc == 1) then
         util_csl_file = name(1:k-1)//".lsj"//".c"
         call openfl(56,util_csl_file,"formatted","UNKNOWN",io_status)
         if (io_status == 1) then
            print *, "Error when opening",util_csl_file
            stop
         end if
      end if
!
!     opening the files  *.lsj.lbl
!
      util_csl_file = name(1:k-1)//".lsj"//".lbl"
      call openfl(57,util_csl_file,"formatted","UNKNOWN",io_status)
      if (io_status == 1) then
         print *, "Error when opening",util_csl_file
         stop
      end if
!GG-2017
      if (unique == 1) then
         util_csl_file = name(1:k-1)//".uni.lsj"//".lbl"
         call openfl(73,util_csl_file,"formatted","UNKNOWN",io_status)
         if (io_status == 1) then
            print *, "Error when opening",util_csl_file
            stop
         end if
         util_csl_file = name(1:k-1)//".uni.lsj"//".sum"
         call openfl(72,util_csl_file,"formatted","UNKNOWN",io_status)
         if (io_status == 1) then
            print *, "Error when opening",util_csl_file
            stop
         end if
      end if
!GG-2017 end
      write(57,"(A53)") " Pos   J   Parity      Energy Total      Comp. of ASF"
!
!     Opening the files   *.lsj.j and
!
      if(ioutj == 1) then
         util_csl_file = name(1:k-1)//".lsj"//".j"
         call openfl(58,util_csl_file,"formatted","UNKNOWN",io_status)
         if (io_status == 1) then
            print *, "Error when opening",util_csl_file
            stop
         end if
         write (58,"(2X,A6,A,F5.1,A,I3,A,I7)" )    &
         name(1:k-1),"  Z = ",nuclear_charge ,"  NEL = ",electron_count,"   NCFG ="! ,asf_set_LS%csf_set_LS%nocsf
      end if
!
!     Opening the files   *.lsj.T and
!
      if(ioutt == 1) then
         util_csl_file = name(1:k-1)//".lsj"//".T"
         open(59,file=util_csl_file,form="unformatted",status="NEW",   &
                                                           iostat=io_status)
         if (io_status /= 0) then
            print *, "Error when opening ",util_csl_file
            stop
         end if
         write (59) "jj2lsj"
         write (59) electron_count, csf_count, orbital_count, nblock
         do jb = 1, nblock
            iatjpdum = two_j(jb) +1
            write (59) jb, ncfinblk(jb), nevinblk(jb), iatjpdum
         end do
      else if(ioutt == 2) then
         util_csl_file = name(1:k-1)//".lsj"//".T"
         open(59,file=util_csl_file,form="unformatted",status="OLD",   &
                                                           iostat=io_status)
         if (io_status /= 0) then
            print *, "Error when opening ",util_csl_file
            stop
         end if
         read (59, iostat=io_status) mixing_file_signature
         if (io_status/=0 .or. mixing_file_signature/="jj2lsj") then
            write (stderr_unit, *) "Not a jj2lsj Transformation File;"
            close(59)
            stop
         end if
         read (59) nelecdum, ncftotdum, nwdum, nblockdum
         if(nelecdum /= electron_count .or. ncftotdum /= csf_count .or. nwdum /= orbital_count   &
                                        .or. nblockdum /= nblock) then
            print*, electron_count, csf_count, orbital_count, nblock
            print*, nelecdum, ncftotdum, nwdum, nblockdum
            print*, "Wrong transformation file *.lsj.T"
            close(59)
            stop
         end if
         do jb = 1, nblockdum
            iatjpdum = two_j(jb) +1
           read (59) jbdum, ncfinblkdum, nevinblkdum, iatjpdum
            if(jbdum /= jb .or. ncfinblkdum /= ncfinblk(jb) .or.       &
               nevinblkdum /= nevinblk(jb) .or.                        &
                                         iatjpdum /= two_j(jb) +1) then
               print*, jb,ncfinblk(jb), nevinblk(jb), two_j(jb) +1
               print*, jbdum, ncfinblkdum, nevinblkdum, iatjpdum
               print*, "Wrong transformation file *.lsj.T"
               close(59)
               stop
            end if
         end do
      end if
      end subroutine inscreen
!
!***********************************************************************
!                                                                      *
      subroutine inscreenlev(record,levels,number_of_levels,fail)
!                                                                      *
!     Attempts to interprete the serial level numbers which are given  *
!     in record. These level numbers can be given in the format:       *
!                1 3 4  7 - 20  48  69 - 85                            *
!     Any order and 'overlapping' intervals are also supported.        *
!     The procedure returns with fail = .true. if the level numbers    *
!     cannot be interpreted properly (fail = .false. otherwise).       *
!                                                                      *
!     The level numbers are returned in the vector                     *
!                                          levels(1:number_of_levels). *
!     The procedure assumes that this vector has a sufficient          *
!     dimension to store all level numbers.                            *
!                                                                      *
!     Calls: idigit.                                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use idigit_i, only: idigit
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      character(len=*), intent(in) :: record
      logical, intent(out)         :: fail
      integer, intent(out)         :: number_of_levels
      integer, intent(inout) :: levels(1:vectors_number)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!GG      logical, dimension(200)   :: low_to
      logical   :: low_to(2000)
      character(len=500)        :: string
      integer                   :: a, i, lower, n
      integer   :: low(2000)
!GG      integer, dimension(200)   :: low
      integer :: run(vectors_number)
!-----------------------------------------------
      fail = .true.;  levels(1:vectors_number) = 0;  number_of_levels = 0
      run(1:vectors_number) = 0;   string = adjustl(record)
      n = 0;   lower = 0;   low_to(:) = .false.
    1 string = adjustl(string)
      if (string(1:1) == " ") then
         if (n == 0  .or.  low_to(n)) then
            return
         else
            goto 10
         end if
      else if (string(1:1) == "-") then
         if (n == 0)  return
         low_to(n)   = .true.;    string(1:1) = " "
         goto 1
      else if (string(1:1) == "0"   .or.   string(1:1) == "1" .or.  &
               string(1:1) == "2"   .or.   string(1:1) == "3" .or.  &
               string(1:1) == "4"   .or.   string(1:1) == "5" .or.  &
               string(1:1) == "6"   .or.   string(1:1) == "7" .or.  &
               string(1:1) == "8"   .or.   string(1:1) == "9") then
         a = idigit(string(1:1));    lower = 10*lower + a
         if (string(2:2) == " "   .or.   string(2:2) == "-") then
            n = n + 1;    low(n) = lower;    lower = 0
         end if
         string(1:1) = " "
         goto 1
      end  if
!
!     Determine no_eigenpairs and max_eigenpair
   10 levels(:) = 0
      do  i = 1,n
         if (low_to(i)) then
            if (low(i) <= low(i+1)) then;   run(low(i):low(i+1)) = 1
            else;                           run(low(i+1):low(i)) = 1
            end if
            if (low_to(i+1)) then;   return
            else;                    cycle
            end if
         else
            run(low(i)) = 1
         end if
      end do
      number_of_levels = 0
      do  i = 1,vectors_number
         if (run(i) == 1) then
            number_of_levels         = number_of_levels + 1
            levels(number_of_levels) = i
         end if
      end do
      fail = .false.
      end subroutine inscreenlev
!
!***********************************************************************
!                                                                      *
      subroutine jj2lsj
      use, intrinsic :: iso_fortran_env, only: real32
!                                                                      *
!     Controls the transformation of atomic states from a jj-          *
!     to a LS-coupled CSF basis.                                       *
!                                                                      *
!     Calls: asf2LS, convrt_double, dallocASFLS, inscreen, ispar,      *
!            itjpo, packlsCSF, prCSFLSall, prCSFjj, prCSFLS, prCSFall, *
!            setLS.                                                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!     Modified by G. Gaigalas,                              May 2021   *
!     Modified by C. Cychen                                     2021   *
!     Modified by G. Gaigalas                                   2022   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: int64
      use eigenstate_expansions,          only: hamiltonian_energy_shift, shifted_eigenvalues, csf_mixing_coefficients
      use orbital_csf_basis,           only: csf_count
      use selected_level_indices,          only: ivec, nvec
      use standard_io_units,        only: stderr_unit
      use numeric_constants,          only: zero, one
      use atomic_calculation_state,           only: nuclear_charge, electron_count
      use transition_basis_blocks,           only: nevinblk, ncfinblk, nblock, two_j
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: itjpo, ispar
      use half_integer_format, only: convrt_double
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!GG NIST
      integer :: i, j, jj, string_l,iblknum,ioutc,ioutj,unique,ioutt
      integer :: level, nocsf_min, lev, string_length
      integer :: nocsf_max, sum_nocsf_min, before_j
!GG NIST
      integer :: ncfmin, ncfmax, ncf_ls_jj_max
      type(ls_ranking), allocatable :: rankings(:)
      integer :: post_lev, post_workers
      real(real64) :: rank_start
      real(real64) :: phase_start, format_start, format_seconds
      real(real64) :: thresh, wa, wb,ssms,g_j,g_jls, sumthrsh
      character(len=4)   :: string_cnum
      character(len=164)  :: string_csf_one
      integer      :: number_of_levels(blocks_number)
      integer :: levels(blocks_number, vectors_number)
      integer     :: max_comp(1:vectors_number)
      integer, pointer     :: leading_ls(:)
      integer, pointer       :: iw(:)
      real(real64), pointer  :: weights(:), weights2(:)
      integer, pointer  :: ithresh(:)
      character(len=164) :: string_csf(1:vectors_number)
      integer(int64) :: loc_i8, ncf_i8
!-----------------------------------------------
      ssms = zero;   g_j = zero;   g_jls = zero;    before_j = 0
      call inscreen(thresh,levels,number_of_levels,ioutc,ioutj,unique, &
                                                                 ioutt)
      allocate(ithresh(csf_count))
      ncf_i8 = csf_count
      do  iblknum = 1, nblock
         if(iblknum == 1) then
            ncfmin = 1
            ncfmax = ncfinblk(iblknum)
         else
            ncfmin = ncfmax + 1
            ncfmax = ncfmin + ncfinblk(iblknum) - 1
         end if
         if(number_of_levels(iblknum) == 0) goto 1
         phase_start=profile_time()
         ithresh = 0
         do  i = ncfmin, ncfmax
            sumthrsh = zero
            do  lev = 1, number_of_levels(iblknum)
                level = levels(iblknum,lev)
                loc_i8 = (level-1)*ncf_i8
                sumthrsh = sumthrsh + dabs(csf_mixing_coefficients(i+loc_i8))
            end do
            if(dabs(sumthrsh) >= dabs(epsnew)) ithresh(i) = 1
         end do
         call profile_emit("screen",phase_start,iblknum,items=count(ithresh(ncfmin:ncfmax)==1))
         phase_start=profile_time()
         write(*,*)"The program generates a list of CSFs in LS-coupling"
         call setls(ithresh,ncfmin,ncfmax)
         call profile_emit("ls_basis",phase_start,iblknum,items=asf_set_ls%csf_set_ls%nocsf)
         phase_start=profile_time()
!
!     output to *.lsj.c
         if(ioutc == 1) call prcsflsall (56,asf_set_ls%csf_set_ls,before_j)
         call profile_emit("csf_write",phase_start,iblknum)
         phase_start=profile_time()
         allocate(asf_set_ls%asf(1:nvec))
         do i = 1, nvec
            allocate(asf_set_ls%asf(i)%eigenvector(1:asf_set_ls%csf_set_ls%nocsf))
         end do
         do lev = 1, number_of_levels(iblknum)
            level = levels(iblknum,lev)
            asf_set_ls%asf(level)%level_no = level
            asf_set_ls%asf(level)%energy = hamiltonian_energy_shift+shifted_eigenvalues(level)
            asf_set_ls%asf(level)%eigenvector = zero
         end do
         asf_set_ls%noasf = number_of_levels(iblknum)
         nocsf_max = asf_set_ls%csf_set_ls%nocsf
!GG NIST
         ncf_ls_jj_max = max(1000,ncfmax-ncfmin+1,asf_set_ls%csf_set_ls%nocsf+1)
         allocate (weights(ncf_ls_jj_max))
         allocate (weights2(ncf_ls_jj_max))
         allocate (iw(ncf_ls_jj_max))
         allocate (leading_ls(ncf_ls_jj_max))
         ncf_i8 = csf_count
         call profile_emit("workspace",phase_start,iblknum)
         do  lev = 1, number_of_levels(iblknum)
            level = levels(iblknum,lev)
            phase_start=profile_time()
            weights = zero;    iw = 0;    wb = zero
            do  i = ncfmin, ncfmax
               if(ithresh(i) == 1) then
                  loc_i8 = (level-1)*ncf_i8
                  wa = csf_mixing_coefficients(i+loc_i8) * csf_mixing_coefficients(i+loc_i8)
                  wb = wb + wa
                  do  j = 1,999
                     if (wa > weights(j)) then
                        weights(j+1:1000) = weights(j:999)
                        weights(j)       = wa
                        iw(j+1:1000)      = iw(j:999)
                        iw(j)            = i - ncfmin + 1
                        exit
                     end if
                  end do
               end if
            end do
            call profile_emit("jj_rank",phase_start,iblknum,level)
            if(lev == 1) then
               print *, " "
               write(*,"(A,A)") " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .",&
               "  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  ."
               write(*, "(A,2X,I4,16X,A,I3)")                            &
               " Under investigation is the block:",iblknum,             &
               " The number of eigenvectors:", number_of_levels(iblknum)
               write(*,"(A,I10,10X,A,I10)")                              &
               " The number of CSF (in jj-coupling):",ncfinblk(iblknum), &
               " The number of CSF (in LS-coupling):",asf_set_ls%csf_set_ls%nocsf
            else
               print *, " "
               write(*,"(A)") " .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  ."
               write(*,"(A)") " The new level is under investigation."
            end if
!
!     perform the transformation
!
            if(lev == 1) then
               write(*,*)                                              &
                "The program generates a list of transformation matrix"
               phase_start=profile_time()
               call asf2ls                                             &
               (iw(1),ithresh,number_of_levels(iblknum),iblknum,levels,&
               ncfmin,ncfmax,ioutt)
               call profile_emit("transform",phase_start,iblknum)
               allocate(rankings(number_of_levels(iblknum)))
               post_workers=post_worker_count(number_of_levels(iblknum))
               phase_start=profile_time()
!$omp parallel do default(none) schedule(static) num_threads(post_workers) if(post_workers>1) &
!$omp shared(rankings,number_of_levels,IBLKNUM,levels,asf_set_LS) private(post_lev,rank_start)
               do post_lev=1,number_of_levels(iblknum)
                  rank_start=profile_time()
                  call rank_ls(asf_set_ls%asf(levels(iblknum,post_lev))%eigenvector,rankings(post_lev))
                  rankings(post_lev)%seconds=profile_time()-rank_start
               end do
!$omp end parallel do
               call profile_emit("ls_rank_wall",phase_start,iblknum,items=post_workers)
               do post_lev=1,number_of_levels(iblknum)
                  call profile_emit("ls_rank_asf",phase_start,iblknum,levels(iblknum,post_lev), &
                       asf_set_ls%csf_set_ls%nocsf,elapsed=rankings(post_lev)%seconds)
               end do
            end if
!
!     output to the screen jj- coupling
            phase_start=profile_time()
            print *, "Weights of major contributors to ASF in jj-coupling:"
            print *, " "
            print *, " Level  J Parity      CSF contributions"
            print *, " "
            if (level > 1   .and.   dabs(wb) > 1.0001_real64) then
               print *, "level, wb = ",level,wb
               stop "JJ2LSJ(): program stop A."
            end if
            nocsf_min = 5
            do j = 1,5
               if(abs(weights(j)) < 0.00001_real64) then
                  nocsf_min = j - 1
                  exit
               end if
            end do
!CYC-2023 Corrections, fix the wrong output information.
!            IF(ISPAR(iw(1)) == -1) THEN
            if(ispar(iw(1)+ncfmin-1) == -1) then
              asf_set_ls%asf(level)%parity = "-"
!            ELSE IF(ISPAR(iw(1)) == 1) THEN
            else if(ispar(iw(1)+ncfmin-1) == 1) then
!CYC-2023 end
              asf_set_ls%asf(level)%parity = "+"
            else
              stop "JJ2LSJ: program stop D."
            end if
            asf_set_ls%asf(level)%totalj = itjpo(iw(1)+ncfmin-1) - 1
            call convrt_double(asf_set_ls%asf(level)%totalj,string_cnum,string_l)
            print '(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7/16x))', &
               ivec(level),string_cnum(1:string_l),               &
              asf_set_ls%asf(level)%parity,(weights(j),iw(j),j=1,nocsf_min)
            print*, "              Total sum over  weight (in jj) is:",wb
            print *, " "
            print *, "Definition of leading CSF:"
            print *, " "
!CGG            call prCSFjj(-1,iw(1))
            call prcsfjj(-1,iw(1)+ncfmin-1,iw(1))
!
!     output to the screen LS- coupling
            print *, " "
            print *, " "
            print *, "Weights of major contributors to ASF in LS-coupling:"
            print *, " "
            print *, " Level  J Parity      CSF contributions"
            print *, " "
            call profile_emit("jj_stdout",phase_start,iblknum,level)
            phase_start=profile_time()
            sum_nocsf_min = 0
            weights = zero;  weights2 = zero;  iw = 0;  wb = zero
            wb=rankings(lev)%norm
            do i=1,asf_set_ls%csf_set_ls%nocsf
               iw(i)=rankings(lev)%order(i)
               weights(i)=asf_set_ls%asf(level)%eigenvector(iw(i))
               weights2(i)=weights(i)*weights(i)
            end do
            call profile_emit("ls_rank_gather",phase_start,iblknum,level,asf_set_ls%csf_set_ls%nocsf)
            phase_start=profile_time()
            if(imincompoff == 1) then
               if (level > 1 .and. (wb*100.0_real32+dabs(mincomp)) <= 100.0_real32) then
                  write (*,"(A,F5.0,A)")   &
                    "Attention!!!  The program not reach the accuracy",mincomp," %"
                  print *, "level, sum of weights in LSJ = ",level,wb
                  print *,""
               end if
            else
               if (level > 1 .and. abs(wb) > 1.001_real64) then
                  print*,"Attention!!!  The sum of weights is bigger then 1.0001"
                  print *, "level, sum of weights in LSJ = ",level,wb
               end if
               if (level > 1 .and.  abs(wb) < 0.998_real64) then
                  print*,"Attention!!!  The sum of weights is little then .998"
                  print *, "level, sum of weights in LSJ = ",level,wb
               end if
            end if
            nocsf_min = nocsf_max;   jj = 0
            do j = 1,nocsf_max
               if(dabs(weights2(j)) < dabs(thresh)) then
                  nocsf_min = j - 1
                  exit
               end if
               jj = jj + 1
               leading_ls(jj) = iw(j)
            end do
            max_comp(level) = iw(1)
            sum_nocsf_min = sum_nocsf_min + nocsf_min
            if (nocsf_min <= 0) then
               stop "JJ2LSJ(): program stop C."
            end if
            call profile_emit("candidate_select",phase_start,iblknum,level,nocsf_min)
            phase_start=profile_time()
            call convrt_double(asf_set_ls%csf_set_ls%csf(iw(1))%totalj,string_cnum,string_l)
            print '(1x,i4,1x,2a4,2x,100(3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7,3x,f8.5," of",i7/16x))', &
               asf_set_ls%asf(level)%level_no,string_cnum(1:string_l), &
                   asf_set_ls%csf_set_ls%csf(iw(1))%parity,                   &
                   (weights2(j),iw(j),j=1,5)
            print*, "              Total sum over  weight (in LSJ) is:",wb
            print *, " "
            print *, "Definition of leading CSF:"
            print *, " "
            if(sum_nocsf_min > 5) sum_nocsf_min = 5
            do  i = 1,csf_count
               do  j = 1, sum_nocsf_min
                  if (i == iw(j)) then
                     call prcsfls (-1,asf_set_ls%csf_set_ls,leading_ls(j))
                     exit
                  end if
               end do
            end do
!
!     output to *.lsj.lbl
            call profile_emit("ls_stdout",phase_start,iblknum,level)
            phase_start=profile_time()
            format_seconds=0.0_real64
            if(before_j == 0) then
               before_j = 1
            else if(lev == 1 .and. before_j == 1) then
               write(57,"(A1)") " "
            end if
            do j = 1,nocsf_min
               format_start=profile_time()
               call packlscsf(asf_set_ls%csf_set_ls,iw(j),string_csf_one)
               format_seconds=format_seconds+(profile_time()-format_start)
               if(j == 1) then
                  write(57,"(I3,1X,A4,5X,A1,8X,F16.9,5X,F7.3,A)")        &
                  asf_set_ls%asf(level)%level_no,string_cnum(1:string_l),&
                  asf_set_ls%csf_set_ls%csf(iw(1))%parity,               &
                  asf_set_ls%asf(level)%energy,wb*100,"%"
                  string_csf(level) = string_csf_one
               end if
               string_length = len_trim(string_csf_one)
               write(57,"(7X,F12.8,3X,F11.8,3X,A)") weights(j),weights2(j),string_csf_one(1:string_length)
            end do
!     output  *.lsj.j and
            call profile_emit("label_format",phase_start,iblknum,level,nocsf_min,elapsed=format_seconds)
            call profile_emit("label_write",phase_start,iblknum,level, &
                              elapsed=profile_time()-phase_start-format_seconds)
            phase_start=profile_time()
            if(ioutj == 1) then
               if(lev == 1) then
                  write (58, "(//A8,I4,2X,A8,I4)" ) "  2*J = ",       &
                  asf_set_ls%asf(level)%totalj,"NUMBER =",number_of_levels(iblknum)
               end if
               string_length = len_trim(string_csf(level))
               write(58,"(3(A8,F15.10))") "Ssms=", ssms, "g_J=",g_j, "g_JLS=", g_jls
               write(58,"(I6,F16.9,2X,A)")                                &
                 max_comp(level),asf_set_ls%asf(level)%energy,string_csf(level)(1:string_length)
               write(58,"(7F11.8)")                                       &
                 (asf_set_ls%asf(level)%eigenvector(i), i=1,asf_set_ls%csf_set_ls%nocsf)
            end if
            call profile_emit("vector_write",phase_start,iblknum,level)
         end do
         phase_start=profile_time()
         call dallocasfls(asf_set_ls)
         deallocate(weights)
         deallocate(weights2)
         deallocate(iw)
         deallocate(leading_ls)
         deallocate(rankings)
         call profile_emit("release",phase_start,iblknum)
    1    continue
      end do
      deallocate(ithresh)
      if(ioutc == 1) call profile_file_size("csf_bytes",56)
      close(56)
!GG-2017
      if (unique == 1) then
         phase_start=profile_time()
         call uniquelsj
         call profile_emit("unique_total",phase_start)
         call profile_file_size("unique_summary_bytes",72)
         call profile_file_size("unique_label_bytes",73)
         close(72)
         close(73)
      end if
!GG-2017
      call profile_file_size("label_bytes",57)
      close(57)
      if(ioutj == 1) write(58,"(A)") "**"
      if(ioutj == 1) call profile_file_size("vector_bytes",58)
      close(58)
      deallocate(jcoup)
      return
      end subroutine jj2lsj
!
!***********************************************************************
!                                                                      *
      subroutine packlscsf(csf_set_ls,csf_index,string_csf)
!                                                                      *
!     Encording all CSF in LS - coupling.                              *
!                                                                      *
!     Rules for encoding                                               *
!        1. All blanks deleted                                         *
!        2. If Qi=1, omit Qi                                           *
!        3. If Qi=1 or Qi>=4l+1, omit ALFAi                            *
!        4. If i=1 or (Qi=4l+2 and i<>m), insert '.'; else _BETAi.     *
!                                                                      *
!     Calls: getchLS, packls_lbl, packls_j.                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
      use packls_i, only: packls
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      type(csf_basis_ls), intent(in) :: csf_set_ls
      integer, intent(in)            :: csf_index
!GG      CHARACTER(LEN=64), INTENT(OUT) :: string_CSF
      character(len=164), intent(out) :: string_csf
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                 :: i, j, counter
      integer, pointer          :: occupation(:)
!GG      integer, dimension(8)                   :: Q
      integer                   :: q(18)
      character(len=4)                        :: ls, xls
      character(len=4),pointer  :: string_ls(:), string_xls(:)
!GG      character(len=3), dimension(15)         :: COUPLE
      character(len=3)         :: couple(35)
!GG      character(len=3), dimension(8)          :: ELC
      character(len=3)          :: elc(18)
      character(len=2)                        :: string
!GG      CHARACTER(LEN=64)                       :: string_CSF_PACK
      character(len=164)                       :: string_csf_pack
      character(len=3), pointer :: string_nl(:)
      character(len=1)       :: l1(0:20)
      data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
               "o","q","r","t","u","v","w","x","y","z"/
!-----------------------------------------------
      counter = 0
      allocate(occupation(csf_set_ls%nwshells))
      allocate(string_ls(csf_set_ls%nwshells))
      allocate(string_xls(csf_set_ls%nwshells))
      do j = csf_set_ls%nwcore+1, csf_set_ls%nwshells
         if (csf_set_ls%csf(csf_index)%occupation(j) > 0) then
            counter = counter +1;        occupation(counter) = j
            call getchls_lbl(csf_index,j,ls,xls)
            string_ls(counter)  = ls;    string_xls(counter) = xls
         end if
      end do
!
      allocate(string_nl(counter))
      do i = 1,counter
         write(string,"(I2)")csf_set_ls%shell(occupation(i))%n
         string_nl(i)(1:2) = string(1:2)
         string_nl(i)(3:3) = l1(csf_set_ls%shell(occupation(i))%l)
      end do
!GG      IF(counter > 8) THEN
      if(counter > 18) then
            print*, counter
            stop "prCSFLS: A"
      end if
      do i = 1,counter
         couple(i)(1:3) = string_ls(i)(1:3)
         if(i > 1) then
            couple(i+counter-1)(1:3) = string_xls(i)(1:3)
         end if
         q(i) = csf_set_ls%csf(csf_index)%occupation(occupation(i))
         elc(i)(1:3) = string_nl(i)(1:3)
      end do
      call packls(counter,elc,q,couple,string_csf_pack)
      string_csf = string_csf_pack
      deallocate(string_nl)
      deallocate(occupation)
      deallocate(string_ls)
      deallocate(string_xls)
      end subroutine packlscsf
!
!***********************************************************************
!                                                                      *
      subroutine prcsfjj(stream,csf_index,csf_number_local)
!                                                                      *
!     Print all information about the single CSF scheme csf_set in     *
!     a nead format on stream.                                         *
!                                                                      *
!     Calls: convrt_double, itjpo, iq, jqs.                            *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis,           only: orbital_count, orbital_principal_n, orbital_kappa
      use csf_pair_angular_state,             only: core_orbital_count => ncore
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: itjpo, iq, jqs
      use half_integer_format, only: convrt_double
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: stream, csf_index, csf_number_local
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                 :: i, orbital_index, counter, text_length
      integer, pointer          :: occupation(:)
      character(len=4), pointer :: string_j(:), string_x(:)
      character(len=4)                        :: string_cnum
      character(len=2)                        :: string
      character(len=4), pointer ::  string_nl(:)
      character(len=1)       :: l1(0:20)
      data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
               "o","q","r","t","u","v","w","x","y","z"/
!-----------------------------------------------
      counter = 0
      allocate(occupation(orbital_count))
      allocate(string_j(orbital_count))
      allocate(string_x(orbital_count))
      do orbital_index = core_orbital_count+1, orbital_count
         if (iq(orbital_index,csf_index) > 0) then
            counter = counter +1
            occupation(counter) = orbital_index
            if (jqs(3,orbital_index,csf_index)-1 == 0) then
               string_j(counter) = "    "
            else
               call convrt_double(jqs(3,orbital_index,csf_index)-1,string_cnum,text_length)
               string_j(counter) = string_cnum(1:text_length)
            end if
            if (jcoup(orbital_index,csf_index) == 0) then
               string_x(counter) = "     "
            else
               call convrt_double(jcoup(orbital_index,csf_index),string_cnum,text_length)
               string_x(counter) =  string_cnum(1:text_length)
            end if
         end if
      end do
!
      allocate(string_nl(counter))
      do i = 1,counter
         write(string,"(I2)")orbital_principal_n(occupation(i))
         string_nl(i)(1:2) = string(1:2)
         orbital_index = ((iabs(orbital_kappa(occupation(i)))*2)-1+    &
             orbital_kappa(occupation(i))/iabs(orbital_kappa(occupation(i))))/2
         string_nl(i)(3:3) = l1(orbital_index)
         if(orbital_kappa(occupation(i)) <= 0) then
         string_nl(i)(4:4) = " "
         else
         string_nl(i)(4:4) = "-"
         end if
      end do
!
      if (stream == -1) then
         write(*,'(4x,i9,")",100(a4,"(",i2,")",2x))') csf_number_local,(string_nl(orbital_index),                   &
            iq(occupation(orbital_index),csf_index),orbital_index=1,counter)
         write(*,'(18x,a4,100(6x,a4))')(string_j(orbital_index),orbital_index=1,counter)
         call convrt_double(itjpo(csf_index)-1,string_cnum,text_length)
         write(*,'(31x,a4,100(6x,a5))')(string_x(orbital_index),orbital_index=2,counter-1),string_cnum(1:text_length)
      end if
      deallocate(occupation)
      deallocate(string_j)
      deallocate(string_x)
      end subroutine prcsfjj
!
!***********************************************************************
!                                                                      *
      subroutine prcsfls(stream,csf_set_ls,csf_index)
!                                                                      *
!     Print all information about the single CSF scheme csf_set in     *
!     a nead format on stream.                                         *
!                                                                      *
!     Calls: convrt_double, getchLS.                                   *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      use half_integer_format, only: convrt_double
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)            :: stream, csf_index
      type(csf_basis_ls), intent(in) :: csf_set_ls
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                 :: i, j, counter, text_length
      integer, pointer          :: occupation(:)
      character(len=4)                        :: ls, xls
      character(len=4),pointer  :: string_ls(:), string_xls(:)
      character(len=4)                        :: string_cnum
      character(len=2)                        :: string
      character(len=3), pointer :: string_nl(:)
      character(len=1)       :: l1(0:20)
      data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
               "o","q","r","t","u","v","w","x","y","z"/
!-----------------------------------------------
      counter = 0
      allocate(occupation(csf_set_ls%nwshells))
      allocate(string_ls(csf_set_ls%nwshells))
      allocate(string_xls(csf_set_ls%nwshells))
      do j = csf_set_ls%nwcore+1, csf_set_ls%nwshells
         if (csf_set_ls%csf(csf_index)%occupation(j) > 0) then
            counter = counter +1;        occupation(counter) = j
            call getchls(csf_index,j,ls,xls)
            string_ls(counter)  = ls;    string_xls(counter) = xls
         end if
      end do
!
      allocate(string_nl(counter))
      do i = 1,counter
         write(string,"(I2)")csf_set_ls%shell(occupation(i))%n
         string_nl(i)(1:2) = string(1:2)
         string_nl(i)(3:3) = l1(csf_set_ls%shell(occupation(i))%l)
      end do
!
      if (stream == -1) then
         write(*,'(i10,")",4x,100(a3,"(",i2,")",2x))')csf_index,(string_nl(j),                                 &
            csf_set_ls%csf(csf_index)%occupation(occupation(j)),j=1,counter)
         call convrt_double(csf_set_ls%csf(csf_index)%totalj,string_cnum,text_length)
         write(*,'(19x,a4,100(5x,a4))')(string_ls(j),j=1,counter),(string_xls(j),j=2,counter),&
              string_cnum(1:text_length)
      else
         write(stream,'(8(1x,a3,"(",i2,")"))')(string_nl(j),                         &
            csf_set_ls%csf(csf_index)%occupation(occupation(j)),j=1,counter)
         call convrt_double(csf_set_ls%csf(csf_index)%totalj,string_cnum,text_length)
         write(stream,'(1x,15(a4))')(string_ls(j),j=1,counter),(string_xls(j),j=2,counter)
      end if
      deallocate(string_nl)
      deallocate(occupation)
      deallocate(string_ls)
      deallocate(string_xls)
      end subroutine prcsfls
!
!***********************************************************************
!                                                                      *
      subroutine prcsflsall(stream,csf_set_ls,before_j)
!                                                                      *
!     Print all information about the CSF scheme csf_set in a nead     *
!     format on stream.                                                *
!                                                                      *
!     Calls: prCSFLS.                                                  *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)            :: stream
      type(csf_basis_ls), intent(in) :: csf_set_ls
      integer, intent(in)            :: before_j
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer                                  :: i, j
      character (len=2)                        :: string
      character (len=3), pointer :: string_nl(:)
      character (len=1)       :: l1(0:20)
      data l1 /"s","p","d","f","g","h","i","k","l","m","n", &
               "o","q","r","t","u","v","w","x","y","z"/
!-----------------------------------------------
      if(before_j == 0) then
         write(stream,*) " "
         allocate(string_nl(csf_set_ls%nwcore))
         do i = 1,csf_set_ls%nwcore
            write(string,"(I2)")csf_set_ls%shell(i)%n
            string_nl(i)(1:2) = string(1:2)
            string_nl(i)(3:3) = l1(csf_set_ls%shell(i)%l)
         end do
         write(stream,'(18(1x,a3))')(string_nl(j),j=1,csf_set_ls%nwcore)
         deallocate(string_nl)
      end if
!
!     output to *.lsj.c
      do i = 1,csf_set_ls%nocsf
         call prcsfls(stream,csf_set_ls,i)
      end do
      write(stream,*) "*"
      end subroutine prcsflsall
!
!***********************************************************************
!                                                                      *
      subroutine setls(truncation_threshold_index,first_csf_index,last_csf_index)
!                                                                      *
!     This subroutine fills up the variable asf_set_LS%csf_set_LS      *
!     with data generated using the one from asf_set%csf_set.          *
!                                                                      *
!     This subroutine contains the following internal subroutines:     *
!          * subroutine            setLS_action                        *
!          * subroutine            setLS_add_quantum_numbers           *
!          * recursive subroutine  setLS_job_count                     *
!          * function              setLS_equivalent_csfs               *
!                                                                      *
!     Calls: itjpo, iq, setLS_equivalent_csfs, setLS_job_count.        *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use atomic_calculation_state,           only: electron_count
      use orbital_csf_basis,           only: orbital_count, orbital_principal_n, orbital_kappa, csf_count
      use csf_pair_angular_state,             only: core_orbital_count => ncore
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: itjpo, iq
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in)  :: truncation_threshold_index(:)
      integer, intent(in)  :: first_csf_index, last_csf_index
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      type(nl), pointer :: shell_temp(:)
      integer, pointer  :: nllsval(:), nljjval(:)
      type(lsj_list)                  :: nonequiv_csfs_jj
!
      integer :: isubc, isubc2, icsf_jj, icsf_jj2, icsf_jj_real, icsf_ls
      integer :: n, action_type
      logical :: new_one, found_parent_minus, found_parent_plus
!
      integer, pointer :: all_occupation(:)
      integer, pointer :: li(:), si(:), l_i(:), s_i(:), w(:), q(:)
      integer                        :: j
!-----------------------------------------------
      asf_set_ls%csf_set_ls%number_of_electrons = electron_count
!
!     1. define  nl, parent
      allocate(shell_temp(orbital_count));  allocate(nllsval(orbital_count));  allocate(nljjval(orbital_count))
      nljjval = 0;  nllsval = 0
!
      asf_set_ls%csf_set_ls%nwshells = 0
      do isubc = 1, orbital_count, 1
         nljjval(isubc) = orbital_principal_n(isubc)*100 +                              &
                  ((iabs(orbital_kappa(isubc))*2)-1+orbital_kappa(isubc)/iabs(orbital_kappa(isubc)))/2
         new_one = .true.
         do isubc2 = 1, isubc-1, 1
            if(nljjval(isubc) == nljjval(isubc2)) new_one = .false.
         end do
         if(new_one) then
            asf_set_ls%csf_set_ls%nwshells = asf_set_ls%csf_set_ls%nwshells + 1
            shell_temp(asf_set_ls%csf_set_ls%nwshells)%n = orbital_principal_n(isubc)
            shell_temp(asf_set_ls%csf_set_ls%nwshells)%l=                &
            ((iabs(orbital_kappa(isubc))*2)-1+orbital_kappa(isubc)/iabs(orbital_kappa(isubc)))/2
            nllsval(asf_set_ls%csf_set_ls%nwshells) = nljjval(isubc)
         end if
      end do
!
      allocate(asf_set_ls%csf_set_ls%shell(asf_set_ls%csf_set_ls%nwshells))
      allocate(asf_set_ls%csf_set_ls%parent(asf_set_ls%csf_set_ls%nwshells))
!
      do isubc=1, asf_set_ls%csf_set_ls%nwshells, 1
         asf_set_ls%csf_set_ls%shell(isubc) = shell_temp(isubc)
      end do
!
      deallocate(shell_temp)
!
!     begin find parent
      do isubc = 1, asf_set_ls%csf_set_ls%nwshells, 1
         found_parent_minus = .false.;  found_parent_plus = .false.
         do isubc2 = 1, orbital_count, 1
            if(nllsval(isubc) == nljjval(isubc2)) then
               if(orbital_kappa(isubc2) > 0) then
                  found_parent_minus = .true.
                  asf_set_ls%csf_set_ls%parent(isubc)%parent_minus = isubc2
               else
                  found_parent_plus = .true.
                  asf_set_ls%csf_set_ls%parent(isubc)%parent_plus  = isubc2
               end if
            end if
         end do
         if(.not.found_parent_plus)  &
            asf_set_ls%csf_set_ls%parent(isubc)%parent_plus = 0
         if(.not.found_parent_minus)  &
            asf_set_ls%csf_set_ls%parent(isubc)%parent_minus = 0
      end do
!
!     end find parent  --------------------
!
!     2. define the number of "core" shells
!        (the LS shell is supposed to be "core" if:
!        1. l=0 and corresponding jj subshell is "core"
!        2. l<>0 l+ and l - "core" subshells            )
      asf_set_ls%csf_set_ls%nwcore=0
      do isubc=1, asf_set_ls%csf_set_ls%nwshells, 1
         if(asf_set_ls%csf_set_ls%parent(isubc)%parent_minus <= core_orbital_count &
            .and.                                                       &
            asf_set_ls%csf_set_ls%parent(isubc)%parent_plus <= core_orbital_count) &
            asf_set_ls%csf_set_ls%nwcore = asf_set_ls%csf_set_ls%nwcore+1
      end do
!
!     3. form the list of "nonequivalent" csfs_jj
!        (i.e. csfs_jj different in J,parity, or
!         some l's occupation numbers Ni = N_(i+) + N_(i-))
      allocate(nonequiv_csfs_jj%items(csf_count))
      nonequiv_csfs_jj%list_size = 0
      do icsf_jj = first_csf_index, last_csf_index
         new_one = .true.
         if(truncation_threshold_index(icsf_jj) == 1) then
            do icsf_jj2 = first_csf_index , icsf_jj - 1
               if(truncation_threshold_index(icsf_jj2) == 1) then
                  if(setls_equivalent_csfs(icsf_jj,icsf_jj2)) then
                     new_one = .false.
                     exit
                  end if
               end if
            end do
            if(new_one) then
               nonequiv_csfs_jj%list_size = nonequiv_csfs_jj%list_size + 1
               nonequiv_csfs_jj%items(nonequiv_csfs_jj%list_size) = icsf_jj
            end if
         end if
      end do
!
!     4. for each nonequivalent csf_jj find all the csfs_LS
!
!         To avoid the dependency on the number of subshells
!         the recursive subroutine is used
      allocate(li(asf_set_ls%csf_set_ls%nwshells))
      allocate(l_i(asf_set_ls%csf_set_ls%nwshells))
      allocate(si(asf_set_ls%csf_set_ls%nwshells))
      allocate(s_i(asf_set_ls%csf_set_ls%nwshells))
      allocate(q(asf_set_ls%csf_set_ls%nwshells))
      allocate(w(asf_set_ls%csf_set_ls%nwshells))
      allocate(all_occupation(asf_set_ls%csf_set_ls%nwshells))
!
!     4.1 - find the number of csfs_LS
      asf_set_ls%csf_set_ls%nocsf = 0
!
      do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
!
!     set the variable for conviniency ...
         icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
         j = itjpo(icsf_jj_real)-1
!
!     define the occupation numbers
!
         do isubc = 1, asf_set_ls%csf_set_ls%nwshells, 1
            all_occupation(isubc) = 0
            do isubc2 = 1, orbital_count, 1
               if(nllsval(isubc) == nljjval(isubc2))       &
               all_occupation(isubc)=all_occupation(isubc)+iq(isubc2,icsf_jj_real)
            end do  !isubc2
!
         end do  !isubc
         n = 0;    action_type = 1
         call setls_job_count(1, n)
         asf_set_ls%csf_set_ls%nocsf=asf_set_ls%csf_set_ls%nocsf + n
      end do
!
!     4.2 - fill up the csf_LS arrays with the corresponding quantum numbers
!     4.2.1  allocate the arrays
      allocate(asf_set_ls%csf_set_ls%csf(asf_set_ls%csf_set_ls%nocsf))
      do n = 1, asf_set_ls%csf_set_ls%nocsf, 1
         allocate(asf_set_ls%csf_set_ls%csf(n)%occupation(asf_set_ls%csf_set_ls%nwshells))
         allocate(asf_set_ls%csf_set_ls%csf(n)%seniority(asf_set_ls%csf_set_ls%nwshells))
         allocate(asf_set_ls%csf_set_ls%csf(n)%w(asf_set_ls%csf_set_ls%nwshells))
         allocate(asf_set_ls%csf_set_ls%csf(n)%shelll(asf_set_ls%csf_set_ls%nwshells))
         allocate(asf_set_ls%csf_set_ls%csf(n)%shells(asf_set_ls%csf_set_ls%nwshells))
         allocate(asf_set_ls%csf_set_ls%csf(n)%shelllx(asf_set_ls%csf_set_ls%nwshells))
         allocate(asf_set_ls%csf_set_ls%csf(n)%shellsx(asf_set_ls%csf_set_ls%nwshells))
      end do
!
!     4.2.2 - fill them with quantum numbers
      icsf_ls = 0
!
      do icsf_jj=1, nonequiv_csfs_jj%list_size, 1
!
!     set the variable for conveniency ...
         icsf_jj_real = nonequiv_csfs_jj%items(icsf_jj)
         j = itjpo(icsf_jj_real)-1
         do isubc = 1, asf_set_ls%csf_set_ls%nwshells, 1
            all_occupation(isubc)=0
            do isubc2 = 1, orbital_count, 1
               if(nllsval(isubc) == nljjval(isubc2))            &
               all_occupation(isubc)=all_occupation(isubc)+iq(isubc2,icsf_jj_real)
            end do  !isubc2
         end do  !isubc
         action_type = 2
         call setls_job_count(1, n)
      end do
!
      deallocate(nonequiv_csfs_jj%items)
      deallocate(nllsval);  deallocate(nljjval)
!
      deallocate(q);   deallocate(w);    deallocate(all_occupation)
      deallocate(li);  deallocate(l_i);  deallocate(si)
      deallocate(s_i)
!
      contains
!
!***********************************************************************
!                                                                      *
      subroutine setls_action(tip,irez)
!                                                                      *
!     The subroutine defines the "action" of subroutine                *
!     setLS_job_count: whether it counts the number                    *
!     of csfs_LS (asf_set_LS%csf_set_LS%novcsf) or fills the           *
!     arrays of wave functions in LS coupling with                     *
!     asf_set_LS%csf_set_LS%csf(...) with the corresponding            *
!     quantum nubers.                                                  *
!                                                                      *
!     Calls: setLS_add_quantum_numbers.                                *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: tip, irez
!-----------------------------------------------
      if (tip == 1) then
         irez = irez + 1
      else
         call setls_add_quantum_numbers()
      end if
      end subroutine setls_action
!
!***********************************************************************
!                                                                      *
      subroutine setls_add_quantum_numbers()
!                                                                      *
!     The subroutine adds quantum numbers stored in temprorary arrays  *
!     Li, Si, L_i, S_i, w, Q to the corresponding arrays               *
!     of asf_set_LS%csf_set_LS%csf() (i.e. to the arrays of            *
!     corresponding quantum numbers of the wave function in LS         *
!     coupling).                                                       *
!                                                                      *
!     Calls: gettermLS, ispar setLS_action, setLS_job_count.           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: ispar
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isubcx
!-----------------------------------------------
      icsf_ls = icsf_ls + 1
!
      if(icsf_ls  >  asf_set_ls%csf_set_ls%nocsf) then
         stop "setLS_add_quantum_numbers(): program stop A."
      end if
!
      asf_set_ls%csf_set_ls%csf(icsf_ls)%totalj = j
!
      if (ispar(icsf_jj_real) == 1) then
         asf_set_ls%csf_set_ls%csf(icsf_ls)%parity = "+"
      else if (ispar(icsf_jj_real) == -1) then
         asf_set_ls%csf_set_ls%csf(icsf_ls)%parity = "-"
      end if
!
      do isubcx = 1, asf_set_ls%csf_set_ls%nwshells, 1
         asf_set_ls%csf_set_ls%csf(icsf_ls)%occupation(isubcx) =        &
                                             all_occupation(isubcx)
         asf_set_ls%csf_set_ls%csf(icsf_ls)%shelll(isubcx)  = li(isubcx)
         asf_set_ls%csf_set_ls%csf(icsf_ls)%shells(isubcx)  = si(isubcx)
         asf_set_ls%csf_set_ls%csf(icsf_ls)%shelllx(isubcx) = l_i(isubcx)
         asf_set_ls%csf_set_ls%csf(icsf_ls)%shellsx(isubcx) = s_i(isubcx)
         asf_set_ls%csf_set_ls%csf(icsf_ls)%w(isubcx)       = w(isubcx)
         asf_set_ls%csf_set_ls%csf(icsf_ls)%seniority(isubcx) =         &
              2*asf_set_ls%csf_set_ls%shell(isubcx)%l + 1 - q(isubcx)
      end do
      end subroutine setls_add_quantum_numbers
!
!***********************************************************************
!                                                                      *
      recursive subroutine setls_job_count(isubc, rez)
!                                                                      *
!     Recursive subroutine for the calculation of the                  *
!     number of csfs_LS and corresponding quantum numbers.             *
!                                                                      *
!     Calls: setLS_action, setLS_job_count, gettermLS, ittk.           *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: Dec 2015   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use ls_jj_subshell_terms, only: subshell_term_ls
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use angular_triangle_and_parity, only: ittk
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: isubc, rez
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      type(subshell_term_ls) :: ls_terms(120)
      integer :: iterm, nr_terms, suma, delta_j
!GG 2015.12.02 NIST
!GG      integer :: number
      integer :: number, numbergg
      integer :: templmax,templmin,tempsmax,tempsmin,temps,templ
      integer :: l_shell, n
!-----------------------------------------------
      if(isubc>(asf_set_ls%csf_set_ls%nwshells) .or. isubc<1) then
         print *, "isubc = ", isubc
         stop "setLS_job_count(): program stop A."
      end if
!
      if(isubc<=asf_set_ls%csf_set_ls%nwshells) then
         if(all_occupation(isubc)==0) then
            if(isubc>1) then
               li(isubc)  = 0;          l_i(isubc) = l_i(isubc-1)
               si(isubc)  = 0;          s_i(isubc) = s_i(isubc-1)
               if(isubc < asf_set_ls%csf_set_ls%nwshells) then
                  call setls_job_count(isubc + 1, rez)
               else
                  if(ittk(s_i(isubc),l_i(isubc),j)==1)          &
                     call setls_action(action_type, rez) !rez=rez+1
               end if
            else
               li(isubc)  = 0;          l_i(isubc) = 0
               si(isubc)  = 0;          s_i(isubc) = 0
!CYC-2023
! Modification for configurations with empty 1s subshell
               call setls_job_count(isubc + 1, rez)
!CYC-2023 end
            end if
         else
            n = all_occupation(isubc)
            l_shell = asf_set_ls%csf_set_ls%shell(isubc)%l
            call gettermls(l_shell,n,ls_terms,number)
            do iterm=1, number, 1
!GG 2015.12.02 NIST
               n = all_occupation(isubc)
               l_shell = asf_set_ls%csf_set_ls%shell(isubc)%l
               call gettermls(l_shell,n,ls_terms,numbergg)
!GG 2015.12.02 NIST
               li(isubc)  = ls_terms(iterm)%subshell_two_l
               si(isubc)  = ls_terms(iterm)%subshell_two_s
               w(isubc)   = ls_terms(iterm)%w
               q(isubc)   = ls_terms(iterm)%subshell_two_quasispin
               if(isubc==1) then
                  l_i(isubc) = ls_terms(iterm)%subshell_two_l
                  s_i(isubc) = ls_terms(iterm)%subshell_two_s
                  if(asf_set_ls%csf_set_ls%nwshells>1) then
                     call setls_job_count(isubc + 1, rez)
                  else
                     if(ittk(s_i(isubc),l_i(isubc),j)==1)          &
                        call setls_action(action_type, rez) !rez=rez+1
                  end if
               else
                  templmax=l_i(isubc-1)+li(isubc)
                  templmin=abs(l_i(isubc-1)-li(isubc))
                  tempsmax=s_i(isubc-1)+si(isubc)
                  tempsmin=abs(s_i(isubc-1)-si(isubc))
                  do templ = templmin, templmax, 2
                     l_i(isubc) = templ
                     do temps = tempsmin, tempsmax, 2
                        s_i(isubc) = temps
                        if(isubc<asf_set_ls%csf_set_ls%nwshells) then
                           call setls_job_count(isubc+1,rez)
                        else
                           if(ittk(s_i(isubc),l_i(isubc),j)==1)        &
                            call setls_action(action_type, rez) ! rez=rez+1
                        end if
                     end do   ! tempS
                  end do      ! tempL
               end if
            end do            ! iterm
         end if
      end if
      end subroutine setls_job_count
!
!***********************************************************************
!                                                                      *
      function setls_equivalent_csfs(ncsf1,second_csf_index)         result(rez)
!                                                                      *
!     This subroutine defines the "equivalency" of two csfs_jj         *
!     in the sence of generation of csfs_LS                            *
!     number of csfs_LS and corresponding quantum numbers.             *
!                                                                      *
!     Calls: itjpo, ispar, iq.                                         *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis,           only: orbital_count
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: ispar
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer :: ncsf1, second_csf_index
      logical :: rez
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: isubc_ls, orbital_index, nls1, nls2
!-----------------------------------------------
      rez= .true.
      if(itjpo(ncsf1) /= itjpo(second_csf_index)) then
         rez = .false.
      else if (ispar(ncsf1) /= ispar(second_csf_index)) then
         rez = .false.
      else
         do isubc_ls = asf_set_ls%csf_set_ls%nwcore + 1,  &
                       asf_set_ls%csf_set_ls%nwshells, 1
            nls1 = 0;   nls2 = 0
            do orbital_index = 1, orbital_count, 1
               if(nllsval(isubc_ls) == nljjval(orbital_index)) then
                   nls1 = nls1 + iq(orbital_index,ncsf1)
                   nls2 = nls2 + iq(orbital_index,second_csf_index)
               end if
            end do
            if(nls1/=nls2) then
               rez=.false.
               return
            end if
         end do
      end if
      end function setls_equivalent_csfs
!
      end subroutine setls
!
!***********************************************************************
!                                                                      *
      recursive function tralsjj(jj_csf_index,ls_csf_index)          result(wa)
!                                                                      *
!     This procedure return the value of the transformation matrix     *
!     from jj- to LS-coupling scheme in the case of any number of      *
!     open shells.                                                     *
!                                                                      *
!     Calls: coefLSjj, iq, jqs, traLSjjmp.                             *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis,           only: orbital_kappa
      use numeric_constants,          only: zero
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use csf_quantum_queries, only: iq, jqs
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: jj_csf_index,ls_csf_index
      real(real64)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: total_number, shell_number, number_minus, number_plus, &
                 jj_minus, n_minus, q_minus, j_1_minus,                 &
                 jj_plus, n_plus, q_plus, j_1_plus, j_1,                &
                 l_shell, n_ls, w_1, q_1, l_1, s_1
!-----------------------------------------------
      shell_number=asf_set_ls%csf_set_ls%nwcore+1
      number_minus=asf_set_ls%csf_set_ls%parent(shell_number)%parent_minus
      number_plus =asf_set_ls%csf_set_ls%parent(shell_number)%parent_plus
      if (number_minus+1 /= number_plus .and.  &
          number_minus*number_plus  /= 0) then
        stop "traLSjj(): program stop A."
      end if
!
      if (number_minus == 0) then
         jj_minus = iabs(orbital_kappa(number_plus))*2 - 3
         n_minus  = 0;    j_1_minus = 0;   q_minus = (jj_minus + 1)/2
      else
         jj_minus  = iabs(orbital_kappa(number_minus))*2 - 1
         n_minus   = iq(number_minus,jj_csf_index)
         q_minus   = (jj_minus +1)/2 - jqs(1,number_minus,jj_csf_index)
         j_1_minus = jqs(3,number_minus,jj_csf_index)-1
         j_1       = jcoup(number_minus,jj_csf_index)
      end if
!
      if (number_plus == 0) then
         jj_plus = iabs(orbital_kappa(number_minus))*2 + 1
         n_plus  = 0;    j_1_plus = 0;   q_plus = (jj_plus + 1)/2
      else
         jj_plus  = iabs(orbital_kappa(number_plus))*2 - 1
         n_plus   = iq(number_plus,jj_csf_index)
         q_plus   = (jj_plus + 1)/2 - jqs(1,number_plus,jj_csf_index)
         j_1_plus = jqs(3,number_plus,jj_csf_index)-1
         j_1       = jcoup(number_plus,jj_csf_index)
      end if
!
      l_shell=asf_set_ls%csf_set_ls%shell(shell_number)%l
      n_ls   =asf_set_ls%csf_set_ls%csf(ls_csf_index)%occupation(shell_number)
      w_1    =asf_set_ls%csf_set_ls%csf(ls_csf_index)%w(shell_number)
      q_1    =2*l_shell+1-                                              &
              asf_set_ls%csf_set_ls%csf(ls_csf_index)%seniority(shell_number)
!
      l_1     = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shelll(shell_number)
      s_1     = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shells(shell_number)
!
      wa = coeflsjj(l_shell,n_ls,w_1,q_1,l_1,s_1,j_1,                   &
                                  jj_minus,n_minus,q_minus,j_1_minus,   &
                                  jj_plus,q_plus,j_1_plus)
      if (abs(wa) > epsnew) then
         total_number =                                                 &
            asf_set_ls%csf_set_ls%nwshells - asf_set_ls%csf_set_ls%nwcore
         if (total_number == 1) then
            wa = wa
         else if (total_number >= 2) then
            do shell_number = asf_set_ls%csf_set_ls%nwcore + 2,         &
                              asf_set_ls%csf_set_ls%nwshells
               if (abs(wa) > epsnew) then
                  wa = wa * tralsjjmp(shell_number,jj_csf_index,ls_csf_index)
               end if
            end do
         else
            wa = zero
         end if
      else
         wa = zero
      end if
      end function tralsjj
!
!***********************************************************************
!                                                                      *
      recursive function tralsjjmp(shell_number,jj_csf_index,ls_csf_index)   result(wa)
!                                                                      *
!     Return the value of main part of the transformation matrix from  *
!     jj- to LS-coupling scheme in the case of any number of open      *
!     shells.                                                          *
!                                                                      *
!     Calls: coefLSjj, iq, ixjtik, jqs, sixj, nine.                    *
!                                                                      *
!     Written by G. Gaigalas,                                          *
!     NIST                                     last update: May 2011   *
!     Modified by G. Gaigalas                               Jan 2025   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use orbital_csf_basis,           only: orbital_kappa
      use numeric_constants,          only: zero, one
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use wigner_6j_triangle_selection, only: ixjtik
      use csf_quantum_queries, only: iq, jqs
      use wigner_9j, only: nine
      use wigner_6j, only: sixj
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: shell_number,jj_csf_index,ls_csf_index
      real(real64)        :: wa
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(real64) :: wa_sum, rac6, rac9
      integer      :: delta_j, number_minus, number_plus,   &
!GG begin  2025
                      number_plus_1,                        &
!GG end    2025
                      jj_minus, n_minus, q_minus,           &
                      jj_plus, n_plus, q_plus,              &
                      j_i_min, j_i_max, j_i_minus,          &
                      j_i_plus, j_i, j_1_i, jp_1_i, j_1_i1, &
                      l_shell, n_ls, w_i, q_i, l_i, s_i,    &
                      l_1_i, s_1_i, l_1_i1, s_1_i1
!-----------------------------------------------
      wa     = zero;      wa_sum = zero
!
      number_minus = asf_set_ls%csf_set_ls%parent(shell_number)%parent_minus
      number_plus  = asf_set_ls%csf_set_ls%parent(shell_number)%parent_plus
      if (number_minus+1 /= number_plus .and.  &
          number_minus*number_plus  /= 0) then
         stop "tranLSjjmp(): program stop A."
      end if
!
      number_plus_1=asf_set_ls%csf_set_ls%parent(shell_number-1)%parent_plus
      if (number_plus_1 == 0) then
         number_plus_1 =                                            &
                 asf_set_ls%csf_set_ls%parent(shell_number-1)%parent_minus
      end if
!
      if (number_minus == 0) then
         jj_minus  = iabs(orbital_kappa(number_plus))*2 - 1
         n_minus   = 0;    j_i_minus = 0;    q_minus = (jj_minus + 1)/2
         jp_1_i    = jcoup(number_plus_1,jj_csf_index)
      else
         jj_minus  = iabs(orbital_kappa(number_minus))*2 - 1
         n_minus   = iq(number_minus,jj_csf_index)
         q_minus   = (jj_minus + 1)/2 - jqs(1,number_minus,jj_csf_index)
         j_i_minus = jqs(3,number_minus,jj_csf_index)-1
         jp_1_i    = jcoup(number_minus,jj_csf_index)
      end if
!
      if (number_plus == 0) then
         jj_plus  = iabs(orbital_kappa(number_minus))*2 - 1
         q_plus   = (jj_plus + 1)/2
         n_plus   = 0;         j_i_plus = 0;         j_1_i   = jp_1_i
      else
         jj_plus  = iabs(orbital_kappa(number_plus))*2 - 1
         n_plus   = iq(number_plus,jj_csf_index)
         q_plus   = (jj_plus + 1)/2 - jqs(1,number_plus,jj_csf_index)
         j_i_plus = jqs(3,number_plus,jj_csf_index)-1
         j_1_i     = jcoup(number_plus,jj_csf_index)
      end if
!
      j_i_min = iabs(j_i_minus - j_i_plus);   j_i_max = j_i_minus + j_i_plus
      !
      l_shell = asf_set_ls%csf_set_ls%shell(shell_number)%l
      n_ls    = asf_set_ls%csf_set_ls%csf(ls_csf_index)%occupation(shell_number)
      w_i     = asf_set_ls%csf_set_ls%csf(ls_csf_index)%w(shell_number)
      q_i     = 2 * l_shell + 1 -                                          &
                asf_set_ls%csf_set_ls%csf(ls_csf_index)%seniority(shell_number)
!GG begin  2025
      if (n_ls == n_minus + n_plus .and. n_ls == 0) then
         wa = one
      else if (n_ls == n_minus + n_plus .and. n_ls /= 0) then
!GG end    2025
         l_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shelll(shell_number)
         s_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shells(shell_number)
         !
         l_1_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shelllx(shell_number)
         s_1_i = asf_set_ls%csf_set_ls%csf(ls_csf_index)%shellsx(shell_number)
         if (shell_number == 2) then
!GG begin  2025
            j_1_i1 = jcoup(number_plus_1,jj_csf_index)
!GG end    2025
            l_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%shelll(shell_number-1)
            s_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%shells(shell_number-1)
         else
            j_1_i1 =jcoup(number_plus_1,jj_csf_index)
            l_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%shelllx(shell_number-1)
            s_1_i1 =asf_set_ls%csf_set_ls%csf(ls_csf_index)%shellsx(shell_number-1)
         end if
!
         do j_i = j_i_min, j_i_max, 2
            delta_j = ixjtik(j_i_minus,j_i_plus,j_i,j_1_i,j_1_i1,jp_1_i)
            if (delta_j /= 0) then
               call nine(l_1_i1,s_1_i1,j_1_i1,l_i,s_i,j_i,l_1_i,s_1_i,j_1_i,  &
                      1,delta_j,rac9)
               if (delta_j /= 0) then
                 call sixj(j_i_minus,j_i_plus,j_i,j_1_i,j_1_i1,jp_1_i,0,rac6)
                 call nine(l_1_i1,s_1_i1,j_1_i1,l_i,s_i,j_i,l_1_i,s_1_i,j_1_i,&
                      0,delta_j,rac9)
                 wa_sum = wa_sum + (j_i + one) * rac6 * rac9 *                &
                      coeflsjj(l_shell,n_ls,w_i,q_i,l_i,s_i,j_i,              &
                                       jj_minus,n_minus,q_minus,j_i_minus,    &
                                       jj_plus,q_plus,j_i_plus)
               end if
            end if
         end do
         wa = wa_sum * sqrt((jp_1_i+one)*(j_1_i1+one)*(l_1_i+one)*(s_1_i+one))
         if (mod(j_i_minus+j_i_plus+j_1_i1+j_1_i,4) /= 0) wa = - wa
      else
         wa = zero
      end if
      end function tralsjjmp
!
!***********************************************************************
!                                                                      *
      subroutine uniquelsj
      use, intrinsic :: iso_fortran_env, only: iostat_end
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use input_error_handling, only: fail_input_io
!                                                                      *
!     Subroutine defines a unique labels for energy levels             *
!                                                                      *
!     Written by  G. Gaigalas                        Vilnius, May 2017 *
!     Written by  G. Gaigalas                        Vilnius, Feb 2024 *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use iso_fortran_env, only: int32
      implicit none (type, external)
      save  ! Preserve the legacy serial controller workspace.
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
!     The maximum number of levels in the list
      integer, parameter :: lev_no  = 10000
!     The maximum number of mixing coefficients of ASF expension
      integer, parameter :: vec_no  = 10000
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      character(len=256), pointer :: coupling(:, :)
      real(real64), pointer       :: comp(:, :), mix(:, :)
      character(len=3), pointer     :: str_no(:)
      character(len=5), pointer     :: str_j(:)
      character(len=1), pointer     :: str_p(:)
      character(len=256), pointer   :: opt_coupling(:)
      character(len=256), pointer   :: tmp_coupling(:)
      real(real64), pointer         :: energy(:), pro(:)
      real(real64), pointer         :: tmp_comp(:), tmp_mix(:)
      integer, pointer           :: mas_max(:), icount(:), iprgg(:)
      character    :: record*53
      character    :: form*11
      character    :: str*7
      character(len=512) :: input_record
      real(real64) :: phase_start
      real(real64) :: max_comp
      integer      :: io_status, ierr, inum, num_lev, i1,i2,i3,i4
      integer      :: ipr, iprf, num_opt, lev_opt, new_j
      integer(int32)    :: last, length
!-----------------------------------------------
      allocate(coupling(vec_no,lev_no))
      allocate(comp(vec_no,lev_no))
      allocate(mix(vec_no,lev_no))
      allocate(str_no(lev_no))
      allocate(str_j(lev_no))
      allocate(str_p(lev_no))
      allocate(tmp_coupling(vec_no))
      allocate(opt_coupling(lev_no))
      allocate(energy(lev_no))
      allocate(pro(lev_no))
      allocate(tmp_comp(vec_no))
      allocate(tmp_mix(vec_no))
      allocate(mas_max(lev_no))
      allocate(icount(lev_no))
      allocate(iprgg(lev_no))
!
      rewind(57,iostat=io_status)
      if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: rewind labels")
      read (57, "(1A53)", iostat=io_status) record
      if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: label header")
      if (len_trim(record) == 0) call fail_input_io(0,57,"uniquelsj: empty label header")
      write(73, "(1A53)") record
!
      new_j = 0
      do while (new_j < 1)
        phase_start=profile_time()
        num_lev = 0
        new_j = 1
        do
          num_lev = num_lev + 1
          if(num_lev > lev_no) then
            print*, "Please extand the arrays. Now we have Lev_No=",   &
                                                                 lev_no
            call fail_input_io(0,57,"uniquelsj: level capacity exceeded")
          end if
          read(57,"(A)",iostat=io_status) input_record
          if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: level record")
          if (len_trim(input_record) < 50) call fail_input_io(0,57,"uniquelsj: short level record")
          read(input_record,"(A3,A5,5X,A1,8X,F16.9,5X,F7.3)",iostat=io_status)str_no(num_lev), &
              str_j(num_lev),str_p(num_lev),energy(num_lev),pro(num_lev)
          if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: level fields")
          if (.not. ieee_is_finite(energy(num_lev)) .or. .not. ieee_is_finite(pro(num_lev))) &
              call fail_input_io(0,57,"uniquelsj: nonfinite level fields")
          inum = 0
          do
            inum = inum + 1
            if(inum > vec_no) then
              print*, "Please extand the arrays. Now we have Vec_No=", &
                                                                vec_no
              call fail_input_io(0,57,"uniquelsj: component capacity exceeded")
              end if
            read(57,"(A)",iostat=io_status) input_record
            if (io_status == iostat_end) then
               mas_max(num_lev) = inum - 1
               goto 1
            end if
            if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: component record")
            str = input_record(1:7)
            if(str /= "       ") then
               backspace(57,iostat=io_status)
               if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: next level backspace")
               mas_max(num_lev) = inum - 1
             exit
            else if (len_trim(input_record) == 0) then
               mas_max(num_lev) = inum - 1
               new_j = 0
               goto 1
            end if
            if (len_trim(input_record) < 37) call fail_input_io(0,57,"uniquelsj: short component record")
            read(input_record,"(A7,F12.8,3X,F11.8,3X,A)",iostat=io_status)str, &
             mix(inum,num_lev),comp(inum,num_lev),coupling(inum,num_lev)
            if (io_status /= 0) call fail_input_io(io_status,57,"uniquelsj: component fields")
            if (.not. ieee_is_finite(mix(inum,num_lev)) .or. .not. ieee_is_finite(comp(inum,num_lev))) &
                call fail_input_io(0,57,"uniquelsj: nonfinite component")
            if (comp(inum,num_lev) < 0.0_real64) call fail_input_io(0,57,"uniquelsj: negative component weight")
          end do
        end do
   1    continue
        do i1 = 1, num_lev
          if (mas_max(i1) < 1) call fail_input_io(0,57,"uniquelsj: level without components")
          if (comp(1,i1) <= 0.0_real64) call fail_input_io(0,57,"uniquelsj: nonpositive leading weight")
        end do
!
        call profile_emit("unique_read",phase_start,items=num_lev)
        phase_start=profile_time()
        num_opt = 0
        icount = 1
        iprgg = 1
        write(72,"(A)")                                                &
                    "          Composition  Serial No.         Coupling"
        write(72,"(A)") "                       of compos."
        write(72,"(A5,A)") " J = ",trim(str_j(num_lev))
        write(72,"(A)")                                                &
                    "--------------------------------------------------"
        do i1 = 1, num_lev
   2      max_comp = 0.0
          lev_opt = 0
          do i2 = 1, num_lev
            if(icount(i2) == 0) cycle
            if(comp(1,i2) > max_comp) then
               lev_opt = i2
               max_comp = comp(1,i2)
            end if
          end do
          if (lev_opt == 0) call fail_input_io(0,57,"uniquelsj: no positive remaining component")
          num_opt = num_opt + 1
          icount(lev_opt) = 0
          ipr = 1
          if(num_opt == 1) then
            ipr = 1
          else
            i3 = 1
            do while (i3 < num_opt)
             if(trim(opt_coupling(i3))==trim(coupling(ipr,lev_opt)))then
               ipr = ipr + 1
               if(ipr > mas_max(lev_opt)) then
                 call fail_input_io(0,57,"uniquelsj: no unused component label")
               end if
               i3 = 1
             else
               i3 = i3 + 1
             end if
            end do
          end if
          iprgg(lev_opt) = ipr + iprgg(lev_opt) - 1
          if(iprgg(lev_opt) > mas_max(lev_opt)                         &
                                        .and. iprgg(lev_opt) /= 1) then
             print*,                                                   &
             "The program is not able perform the identification for", &
             " level = ",lev_opt
             call fail_input_io(0,57,"uniquelsj: no unique identification")
          end if
          if(ipr == 1) then
            opt_coupling(num_opt) = coupling(ipr,lev_opt)
            iprf = iprgg(lev_opt)
            write(72,"(A,I4,2X,F12.9,I5,3X,A,A,I4)")"Pos",lev_opt,     &
            comp(ipr,lev_opt),iprf,trim(opt_coupling(num_opt))
            if(iprf > 1) then
              tmp_mix = 0
              tmp_comp = 0
              tmp_coupling = ""
              i4 = mas_max(lev_opt)
              tmp_mix(2:iprf)          = mix(i4-iprf+2:i4,lev_opt)
              tmp_mix(iprf+1:i4)       = mix(2:i4-iprf+1,lev_opt)
              mix(2:i4,lev_opt)        = tmp_mix(2:i4)
!
              tmp_comp(2:iprf)         = comp(i4-iprf+2:i4,lev_opt)
              tmp_comp(iprf+1:i4)      = comp(2:i4-iprf+1,lev_opt)
              comp(2:i4,lev_opt)       = tmp_comp(2:i4)
!
              tmp_coupling(2:iprf)    = coupling(i4-iprf+2:i4,lev_opt)
              tmp_coupling(iprf+1:i4) = coupling(2:i4-iprf+1,lev_opt)
              coupling(2:i4,lev_opt)  = tmp_coupling(2:i4)
            end if
          end if
          if(ipr > 1) then
            tmp_mix = 0
            tmp_comp = 0
            tmp_coupling = ""
            i4 = mas_max(lev_opt)
            tmp_mix(1:1) = mix(ipr:ipr,lev_opt)
            tmp_mix(i4-ipr+2:i4)      = mix(1:ipr-1,lev_opt)
            tmp_mix(2:i4-ipr+1)       = mix(ipr+1:i4,lev_opt)
            mix(1:i4,lev_opt)         = tmp_mix(1:i4)
!
            tmp_comp(1:1)             = comp(ipr:ipr,lev_opt)
            tmp_comp(i4-ipr+2:i4)     = comp(1:ipr-1,lev_opt)
            tmp_comp(2:i4-ipr+1)      = comp(ipr+1:i4,lev_opt)
            comp(1:i4,lev_opt)        = tmp_comp(1:i4)
!
            tmp_coupling(1:1)         = coupling(ipr:ipr,lev_opt)
            tmp_coupling(i4-ipr+2:i4) = coupling(1:ipr-1,lev_opt)
            tmp_coupling(2:i4-ipr+1)  = coupling(ipr+1:i4,lev_opt)
            coupling(1:i4,lev_opt)    = tmp_coupling(1:i4)
!
            icount(lev_opt) = 1
            num_opt = num_opt - 1
            goto 2
          end if
        end do
        write(72,"(A)")                                                &
                    "--------------------------------------------------"
        write(72,"(A)")""
!
        call profile_emit("unique_assign_summary",phase_start,items=num_lev)
        phase_start=profile_time()
        do i1 = 1, num_lev
          write(73,"(A3,A5,5X,A1,8X,F16.9,5X,F7.3,A)")                 &
          str_no(i1),str_j(i1),str_p(i1),energy(i1),pro(i1),"%"
          do i2 = 1, mas_max(i1)
            write(73,"(A7,F12.8,3X,F11.8,3X,A)")                       &
            str,mix(i2,i1),comp(i2,i1),trim(coupling(i2,i1))
          end do
        end do
        if(new_j < 1) write(73,"(A)") " "
        call profile_emit("unique_write",phase_start,items=num_lev)
!
      end do
      deallocate(coupling,comp,mix,str_no,str_j,str_p,tmp_coupling,opt_coupling, &
                 energy,pro,tmp_comp,tmp_mix,mas_max,icount,iprgg)
      end subroutine uniquelsj
!
end module jj2lsj_code
