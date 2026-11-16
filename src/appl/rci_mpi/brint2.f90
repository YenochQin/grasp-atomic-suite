!***********************************************************************
!                                                                      *
      subroutine brint2 (ia,ib,ic,id,k,tegral)
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
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      integer :: key, index, loc, ju, jl, jm
!-----------------------------------------------

      key = orbital_count + 1

      index = ((ia*key+ib)*key+ic)*key+id
!
      jl = kstartbreit2(k)
      ju = kstartbreit2(k+1) - 1

      if (index.lt.indtp2(jl).or.index.gt.indtp2(ju)) then
        write(*,*) 'Something wrong in brint2'
        stop
      endif
!
!   The index is within the range of the indices stored; search
!   for it in the list of indices
!
    1 if (ju-jl .gt. 1) then
        jm = (ju+jl)/2
        if (indtp2(jm) .gt. index) then
          ju = jm
        else
          jl = jm
        endif
        goto 1
      endif
!
!   The range is bracketed to the extent possible
!
      if (index .eq. indtp2(ju)) then
        loc = ju
      elseif (index .eq. indtp2(jl)) then
        loc = jl
      else
        write(*,*) k,ia,ib,ic,id,index
        write(*,*) 'Brint2 Integral not found'
        stop
      endif
!
!   Return the value of the integral
!   from storage

      tegral = valtp2(loc)
!
      return
      end subroutine brint2
