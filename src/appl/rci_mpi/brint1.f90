!***********************************************************************
!                                                                      *
      subroutine brint1 (ia,ib,ic,id,k,tegral)
!                                                                      *
!   Returns integrals for the transverse photon interaction.           *
!                                                                      *
!   Written by Per Jonsson                   Octaober 2014             *
!                                                                      *
!***********************************************************************
!...Translated by Charlotte Froese Fischer
!                       Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      use vast_kind_param, only: double
      use bilst_c
      use orb_c,           only: orbital_count
      use kkstartbreit_c
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      integer, intent(in) :: ia, ib, ic, id, k
      real(double), intent(out) :: tegral
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      integer ::  index, loc, ju, jl, jm, key
!-----------------------------------------------
!
      key = orbital_count + 1
      index = ((ia*key+ib)*key+ic)*key+id
!
      jl = kstartbreit1(k)
      ju = kstartbreit1(k+1) - 1
!
      if (index < indtp1(jl).or.index > indtp1(ju)) then
        write(*,*) 'Something wrong in brint1'
        stop
      endif
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 if (ju-jl .gt. 1) then
         jm = (ju+jl)/2
         if (indtp1(jm) > index) then
            ju = jm
         else
            jl = jm
         endif
         goto 1
      endif
!
!   The range is bracketed to the extent possible
!
      if (index .eq. indtp1(ju)) then
         loc = ju
      elseif (index .eq. indtp1(jl)) then
         loc = jl
      else
         write(*,*) k,ia,ib,ic,id,index
         write(*,*) 'Brint1 Integral not found'
         stop
      endif
!
!   Return the value of the integral
!   from storage

      tegral = valtp1(loc)
!
      return
      end subroutine brint1
