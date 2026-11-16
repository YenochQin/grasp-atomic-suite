!***********************************************************************
!                                                                      *
       subroutine rinti_nms (j,k,valnmsk123,valnmsk1)
!                                                                      *
!   The value of this  function is the one-electron integral I (J,K)   *
!   for  orbitals  J, K. The analytical expression for this quantity   *
!   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
!   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
!                                                                      *
!   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
!                                                                      *
!   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
!                                                                      *
!   Modified by C. Naz\'e  Oct. 2011                                   *
!***********************************************************************
!...Translated by Gediminas Gaigalas 11/18/19
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param,  only: double
      use parameter_def,    only: radial_workspace_capacity
      use def_c
      use grid_c
      use wave_c
      use orb_c
      use tatb_c
      use debug_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      use dpbdt_i
      use quad_i
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      real(double), intent(out) :: valnmsk123, valnmsk1
      integer, intent(in) :: j
      integer, intent(in) :: k
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      real(double), dimension(radial_workspace_capacity) :: pa, pb, qa, qb
      real(double) :: a1, a2, piece1, piece2, piece3, piece4
      integer :: l, kap1, kap2, i, l1, l2, j1, j2, l2_tilda
!-----------------------------------------------
!
!
!   Set KAP1 and KAP2
!
!C*      RINTI_NMS = 0.0D 00
      valnmsk1   = 0.0d00
      valnmsk123 = 0.0d00
      kap1 = orbital_kappa(j)
      kap2 = orbital_kappa(k)
!
!   Stop if orbitals J and K have different kappa values
!
      if (kap1 .ne. kap2) then
         write (*,300) orbital_principal_n(j),orbital_symmetry_labels(j),orbital_principal_n(k),orbital_symmetry_labels(k)
         stop
      endif
!
!   Determine the l quantum numbers
!
      if (kap1 .gt. 0) then
         l1 =  kap1
      else
         l1 = -kap1-1
      endif
!
      if (kap2 .gt. 0) then
         l2 =  kap2
      else
         l2 = -kap2-1
      endif
      if (l1 .ne. l2) return
!
!   Determine the j quantum numbers and l_2 tilda
!
      j1 = 2*iabs(orbital_kappa(j))-1
      j2 = 2*iabs(orbital_kappa(k))-1
      l2_tilda = j2 - l2
!
      mtp = max (orbital_last_grid_index(j),orbital_last_grid_index(k))
!
!   Kinetic energy contribution
!
!   Piece involving derivatives
!
      call dpbdt (j)
      pa(1) = 0.0d00
      qa(1) = 0.0d00
      do 11 i = 2,mtp
        pa(i) = ta(i)
        qa(i) = tb(i)
   11 continue
      call dpbdt (k)
      pb(1) = 0.0d00
      qb(1) = 0.0d00
      do 12 i = 2,mtp
        pb(i) = ta(i)
        qb(i) = tb(i)
   12 continue
      ta(1) = 0.0d00
      do 1 i = 2,mtp
         ta(i) = (pa(i)*pb(i)+qa(i)*qb(i))/radial_grid_derivative(i)
    1 continue
      call quad (piece1)
      piece1 = piece1/(transformed_grid_step*transformed_grid_step)
!
      a1 = dble(l2*(l2+1))
      a2 = dble(l2_tilda*(l2_tilda+1))
      ta(1) = 0.0d00
      do 2 i = 2,mtp
        ta(i) = radial_grid_derivative(i)                                                  &
              *(a1*radial_large_components(i,j)*radial_large_components(i,k)+a2*radial_small_components(i,j)* &
                 & radial_small_components(i,k))                 &
              /(radial_grid(i)*radial_grid(i))
    2 continue
      call quad (piece2)
!
!   Pieces not involving derivatives
!
      ta(1) = 0.d0
      do 3 i = 2,mtp
         ta(i) = (radial_small_components(i,j)*pb(i)+radial_small_components(i,k)*pa(i))/radial_grid(i)
    3 continue
      call quad (piece3)
      piece3 = -2.0*nuclear_charge*piece3/(speed_of_light_au*transformed_grid_step)
!
      if (kap2 .ne. 1) then
        ta(1) = 0.0d00
        do 4 i = 2,mtp
           ta(i) = radial_grid_derivative(i)                                               &
                 * (radial_small_components(i,j)*radial_large_components(i,k)+radial_large_components(i,j)* &
                    & radial_small_components(i,k))                   &
                 /(radial_grid(i)*radial_grid(i))
    4   continue
        call quad (piece4)
        piece4 = -piece4*nuclear_charge*(dble (kap2)-1.0)/speed_of_light_au
      else
        piece4 = 0.0d00
      end if
      valnmsk1 = 0.5*(piece1+piece2)
      valnmsk123 = 0.5*(piece1+piece2+piece3+piece4)
!
!   Debug printout
!
!C*      LDBPR(5)=.TRUE.
      if (ldbpr(5))                                                    &
         write (99,302) orbital_principal_n(j),orbital_symmetry_labels(j),orbital_principal_n(k), &
            & orbital_symmetry_labels(k),valnmsk123
!
!C*      RETURN
!
  300 format ('RINTI_NMS: Attempt to calculate I(',                    &
               1i2,1a2,',',1i2,1a2,')')
  302 format (/' K (',1i2,1a2,',',1i2,1a2,') = ',1pd19.12)
!
      return
      end subroutine rinti_nms
