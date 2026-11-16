!***********************************************************************
!                                                                      *
      subroutine setlagmpi(eol)
!                                                                      *
!   Sets up the data structure  pertaining to the Lagrange multipli-   *
!   ers  on the first entry;  on subsequent calls it  determines new   *
!   estimates for the multipliers.                                     *
!                                                                      *
!   Call(s) to: [LIB92]: ALLOC, QUAD, RINTI,                           *
!               [RSCF92]: DACON, SETCOF, XPOT,  YPOT.                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford    Last update: 18 Dec 1992   *
!   MPI version by Xinghong He              Last update: 03 Aug 1998   *
!                                                                      *
!***********************************************************************
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use, intrinsic :: iso_fortran_env, only: real64
      use parameter_def,    only: orbital_index_radix, max_radial_points
      use memory_man
      use orba_c
      use core_c
      use def_c
      use fixd_c
      use grid_c, only: radial_point_count, radial_derivative_over_radius
      use lagr_c
      use orb_c
      use mpi_c
      use pote_c, only: yp, xp, xq
      use scf_c
      use tatb_c, only: mtp, ta
      use wave_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use setcof_i
      use ypot_i
      use xpot_i
      use dacon_i
      use quad_i
      use rinti_i
      implicit none (type, external)
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      logical  :: eol
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      real(real64), parameter :: p001 = 1.0e-1_real64
      integer, parameter :: key = orbital_index_radix
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer :: itwice, liraw, li, lip1, nakli, ljraw, lj, ieccli, l1, l2, &
         jlast, mlast, m, j, i
      real(real64) :: ypj(max_radial_points), ypm(max_radial_points), xpj(max_radial_points), xpm(max_radial_points), &
         xqj(max_radial_points), xqm(max_radial_points)
      real(real64) :: eps,ucfj,ucfm,result,rijm,qdif,obqdif,obqsum,tmp
      logical :: first, fixli, fixlj, fulli, fullj
!-----------------------------------------------
!
      data first/ .true./
!
!-----------------------------------------------------------------------

      if (first) then

!=======================================================================
!   Determine the total number of Lagrange multipliers and store
!   their indeces in IECC(1:NEC). Memories are allocated for IECC
!   and ECV. The only outputs are NEC and IECC(). This is done only
!   for the first time this routine is called.
!
!   Implementation changes:
!   The "DO LIraw ..." is executed twice (if NEC obtained in the
!   first time > 0) so the memory allocation is done only once.
!
!   This part is not distributed.
!=======================================================================

         eps = accy*0.01_real64      ! criterion to see if an orb is occupied
         do itwice = 1, 2
            nec = 0
!            IF (ITWICE /= 2) THEN
!               DO LIRAW = 1, NW - 1
!                  LI = IORDER(LIRAW)
!                  LIP1 = MAX(NCORE,LIRAW) + 1
!                  NAKLI = NAK(LI)
!                  FIXLI = LFIX(LI)
!                  FULLI = ABS(UCF(LI)-DBLE(NKJ(LI)+1)) < EPS
!                  DO LJRAW = LIP1, NW
!                     LJ = IORDER(LJRAW)
!                     FIXLJ = LFIX(LJ)
!                     FULLJ = ABS(UCF(LJ)-DBLE(NKJ(LJ)+1)) < EPS
!                     IF (.NOT.(NAK(LJ)==NAKLI .AND. .NOT.(FIXLI .AND. FIXLJ)&
!                         .AND. .NOT.(FULLI .AND. FULLJ))) CYCLE
!                     NEC = NEC + 1
!                     CYCLE
!                  !*** Encode index at 2nd round ***
!                  END DO
!               END DO
!            ELSE
!               DO LIRAW = 1, NW - 1
!                  LI = IORDER(LIRAW)
!                  LIP1 = MAX(NCORE,LIRAW) + 1
!                  NAKLI = NAK(LI)
!                  FIXLI = LFIX(LI)
!                  FULLI = ABS(UCF(LI)-DBLE(NKJ(LI)+1)) < EPS
!                  !*** Encode index at 2nd round ***
!                  DO LJRAW = LIP1, NW
!                     LJ = IORDER(LJRAW)
!                     FIXLJ = LFIX(LJ)
!                     FULLJ = ABS(UCF(LJ)-DBLE(NKJ(LJ)+1)) < EPS
!                     IF (.NOT.(NAK(LJ)==NAKLI .AND. .NOT.(FIXLI .AND. FIXLJ)&
!                         .AND. .NOT.(FULLI .AND. FULLJ))) CYCLE
!                     NEC = NEC + 1
!                  !*** Encode index at 2nd round ***
!                     IECC(NEC) = LI + KEY*LJ
!                  END DO
!               END DO
!            ENDIF
         do liraw = 1, orbital_count - 1
            li = iorder(liraw)
            lip1 = max (closed_core_subshell_count, liraw) + 1
            nakli = orbital_kappa(li)
            fixli = is_orbital_fixed(li)
            fulli = abs ( ucf(li)-dble (orbital_two_j(li)+1) ) .lt. eps
            do ljraw = lip1, orbital_count
               lj = iorder(ljraw)
               fixlj = is_orbital_fixed(lj)
               fullj = abs ( ucf(lj)-dble (orbital_two_j(lj)+1) ) .lt. eps
               if ( (orbital_kappa(lj) .eq. nakli) .and.            &
                     (.not. (fixli .and. fixlj)) .and.    &
                     (.not. (fulli .and. fullj)) ) then
                  nec = nec + 1
                  !*** Encode index at 2nd round ***
                  if (itwice == 2) iecc(nec) = li + key * lj
               endif
            enddo
         enddo


            if (itwice==1 .and. nec>0) then
               call alloc (ecv, nec, 'ECV', 'SETLAGmpi')
               call alloc (iecc, nec, 'IECC', 'SETLAGmpi')
            else
               exit
            endif
         end do                                  !itwice

!=======================================================================
!   Print information about Lagrange multipliers
!=======================================================================

         if (mpi_rank == 0) then
            if (nec == 0) then
               write (*, 302)
            else
               write (*, 304)
               do li = 1, nec
               !*** Decode index ***
                  ieccli = iecc(li)
                  l1 = ieccli/key
                  l2 = ieccli - key*l1
                  write (*, 305) orbital_principal_n(l2), orbital_symmetry_labels(l2), orbital_principal_n(l1), &
                     & orbital_symmetry_labels(l1)
               end do
            endif
         endif
         first = .false.
      endif

!FF+GG  12/07/05
!     Lagrange multipliers need to be computed also on the first call
!     RETURN

      if (nec == 0) return
      if (mpi_rank == 0) write (*, 306)
      jlast = 0
      mlast = 0

      do li = 1, nec
         !*** Decode index ***
         ieccli = iecc(li)
         m = ieccli/key
         j = ieccli - key*m
!
         if (j /= jlast) then
            ucfj = ucf(j)
            call setcof (eol, j)
            call ypot (j)
            call xpot (j)
            call dacon
            ypj(:radial_point_count) = yp(:radial_point_count)
            xpj(:radial_point_count) = xp(:radial_point_count)
            xqj(:radial_point_count) = xq(:radial_point_count)
            jlast = j
         endif
!
         if (m /= mlast) then
            ucfm = ucf(m)
            call setcof (eol, m)
            call ypot (m)
            call xpot (m)
            call dacon
            ypm(:radial_point_count) = yp(:radial_point_count)
            xpm(:radial_point_count) = xp(:radial_point_count)
            xqm(:radial_point_count) = xq(:radial_point_count)
            mlast = m
         endif
!
         mtp = max(orbital_last_grid_index(j),orbital_last_grid_index(m))
!
         if (is_orbital_fixed(m)) then
            ta(1) = 0._real64
            do i = 2, mtp
               ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                  & -radial_small_components(i,m)*xpj(i))*speed_of_light_au+(radial_large_components(i,m)*radial_large_components(i&
                  ,j)+radial_small_components(i,m)*radial_small_components(i,j))*ypj(i))
            end do

            call quad (result)
            rijm = rinti(m,j,1)
            ecv(li) = (result - rijm / mpi_process_count)*ucfj
! start dbg
!           WRITE (81,*)'1, RESULT, RIJM, UCFJ, ECV, TA' ! dbg
!           WRITE (81,*)RESULT, RIJM, UCFJ, ECV ! dbg
!           DO i = 1, MTP ! dbg
!              WRITE (81,*) i, TA(i), r(i), rp(i) ! dbg
!           ENDDO ! dbg
! end dbg

         else if (is_orbital_fixed(j)) then
            ta(1) = 0._real64
            do i = 2, mtp
               ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,j)*xqm(i) &
                  & -radial_small_components(i,j)*xpm(i))*speed_of_light_au+(radial_large_components(i,j)*radial_large_components(i&
                  ,m)+radial_small_components(i,j)*radial_small_components(i,m))*ypm(i))
!GG             write(222,*)"XQM(I)",XQM(I)
            end do

!start dbg
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i)
!               WRITE (83,*) i, r(i), rp(i), rpor(i)
!               WRITE (84,*) i, pf(i,j), qf(i,j)
!               write(85,*) i,ypm(i)
!              write(86,*)i,xpm(i),xqm(i)
!           ENDDO
! end dbg
            call quad (result)

            rijm = rinti(j,m,1)                  !/ nprocs
            ecv(li) = (result - rijm / mpi_process_count)*ucfm
!start dbg
!           WRITE (81,*)'2, RESULT, RIJM, UCFM, ECV, TA'
!           WRITE (81,*)RESULT, RIJM, UCFJ, ECV, r(i), rp(i)
!end dbg


         else
            qdif = abs((ucfj - ucfm)/max(ucfj,ucfm))
            if (qdif > p001) then
               obqdif = 1._real64/ucfj - 1._real64/ucfm
               ta(1) = 0._real64
               do i = 2, mtp
                  ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                     & -radial_small_components(i,m)*xpj(i)-radial_large_components(i,j)*xqm(i&
                     )+radial_small_components(i,j)*xpm(i))*speed_of_light_au+(ypj(i)-ypm(i))* &
                        & (radial_large_components(i,m)*radial_large_components(i,j)+radial_small_components(i,&
                     m)*radial_small_components(i,j)))
               end do

               call quad (result)
               ecv(li) = result/obqdif
!start dbg
!           WRITE (81,*)'3, RESULT, OBQDIF, ECV, TA'
!           WRITE (81,*)RESULT, OBQDIF, ECV
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i), r(i), rp(i)
!           ENDDO
!end dbg


            else
               obqsum = 1._real64/ucfj + 1._real64/ucfm
               ta(1) = 0._real64
               do i = 2, mtp
                  ta(i) = radial_derivative_over_radius(i)*((radial_large_components(i,m)*xqj(i) &
                     & -radial_small_components(i,m)*xpj(i)+radial_large_components(i,j)*xqm(i&
                     )-radial_small_components(i,j)*xpm(i))*speed_of_light_au+(ypj(i)+ypm(i))* &
                        & (radial_large_components(i,m)*radial_large_components(i,j)+radial_small_components(i,&
                     m)*radial_small_components(i,j)))
               end do

               call quad (result)
               rijm = rinti(m,j,1)               !/ nprocs
               ecv(li) = (result - 2._real64*rijm / mpi_process_count)/obqsum
!start dbg
!           WRITE (81,*)'4, RESULT, RIUJM, OBQSUM, ECV, TA'
!           WRITE (81,*)RESULT, RIUJM, OBQSUM, ECV
!           DO i = 1, MTP
!              WRITE (81,*) i, TA(i), r(i), rp(i)
!           ENDDO
!end dbg

            endif
         endif
!=======================================================================
!  Collect contributions from all nodes.
!  Another alternative is to modify mcpmpi to let every node
!  have the same set of "LAB" and "NCONTR". So an easier aggregation
!  can be done for arrays DA, XA, YA in setcof. And so the computation
!  of YPOT, XPOT, DACON can be distributed.
!=======================================================================
         call mpi_allreduce (ecv(li),tmp,1,mpi_double_precision, &
                          mpi_sum, mpi_comm_world, mpi_error_code)
         ecv(li) = tmp

         if(mpi_rank == 0)write (*, 307) orbital_principal_n(j), orbital_symmetry_labels(j), orbital_principal_n(m), &
            & orbital_symmetry_labels(m), ecv(li)

      end do

!db      close(81)
!db      close(82)


  302 format(/,'Lagrange multipliers are not required')
  304 format(/,'Include Lagrange multipliers between:'/)
  305 format(13x,2(2x,1i2,1a2))
  306 format(/,'Lagrange multipliers:'/)
  307 format(13x,2(2x,1i2,1a2),2x,1p,d16.9)

      return
      end subroutine setlagmpi
