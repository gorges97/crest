!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Philipp Pracht
!
! crest is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! crest is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with crest.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================!

module msmod

      use iso_fortran_env, wp => real64
      implicit none

      
      !-- storage for a single mol
      type msmol
          integer :: nat
          integer,allocatable :: at(:)
          real(wp),allocatable :: xyz(:,:)
          real(wp) :: etot
          integer :: chrg
          integer :: nfrag
          integer,allocatable :: frag(:)
          integer :: gen
          contains
              procedure :: newmol => msmod_newmol
              procedure :: dist => msmoldist
      end type msmol

      public :: msmol

      !-- isomer list
      type ilist
          integer :: nmol
          type(msmol),allocatable :: mol(:)
          logical,allocatable :: new(:)
          contains
              procedure :: dealloc => deallocate_ilist
              procedure :: append => ilist_append
      end type ilist

      !-- fragmentized structure list
      type flist
          integer :: nmol
          type(msmol),allocatable :: mol(:)
          logical,allocatable :: new(:)
          contains
              procedure :: dealloc => deallocate_flist
      end type flist


      !-- global data for msreact
      type msobj

          real(wp) :: ethr    = 1.0_wp  !cov. struc. chekc. thr.
          real(wp) :: wbormsd = 0.5_wp  !wbo rmsd comp. thr.
          
          real(wp) :: e_reac = 0.0_wp ! sum of all reaction (or activation ?) energies "used" to get to this point
          
          real(wp) :: T  = 3000.0_wp    !start temperature
          real(wp) :: fc = 0.05_wp      !start fc
          real(wp) :: cdist = 1.5_wp    !constraing distance scaling factor of rcov
          integer  :: maxc = 30         !max optimization cycle

          type(ilist) :: il   
          type(flist) :: fl
          
      end type msobj

      public :: msobj

      
contains

subroutine deallocate_ilist(self)
     class(ilist) :: self
     if(allocated(self%mol)) deallocate(self%mol)
     if(allocated(self%new)) deallocate(self%new)
     return
end subroutine deallocate_ilist

subroutine deallocate_flist(self)
     class(flist) :: self
     if(allocated(self%mol)) deallocate(self%mol)
     if(allocated(self%new)) deallocate(self%new)
     return
end subroutine deallocate_flist

subroutine msmod_newmol(self,nat,xyz,at,etot,chrg,gen)
     implicit none
     class(msmol) :: self
     integer :: nat
     real(wp) :: xyz(3,nat)
     integer :: at(nat)
     real(wp) :: etot
     integer :: chrg
     integer :: gen
     self%nat = nat
     allocate(self%xyz(3,nat))
     self%xyz = xyz
     allocate(self%at(nat))
     self%at = at
     self%etot = etot
     self%chrg = chrg
     self%gen = gen
     return    
end subroutine msmod_newmol

!calculate the distance between two atoms i and j
function msmoldist(self,i,j) result(dist)
     implicit none
     class(msmol) :: self
     integer :: i,j
     real(wp) :: dist
     dist = 0.0_wp
     if(allocated(self%xyz))then
      dist=(self%xyz(1,i)-self%xyz(1,j))**2 + &
     &     (self%xyz(2,i)-self%xyz(2,j))**2 + &
     &     (self%xyz(3,i)-self%xyz(3,j))**2
      dist = sqrt(dist)
     endif
     return
end function msmoldist



subroutine ilist_append(self,nat,at,xyz,etot,chrg,gen)
     implicit none
     class(ilist) :: self
     type(msmol) :: mol
     integer :: nat
     real(wp) :: xyz(3,nat)
     integer :: at(nat)
     real(wp) :: etot
     integer :: chrg
     integer :: gen

     type(msmol),allocatable :: dummy(:)
     logical,allocatable     :: btmp(:)
     integer :: n
     call mol%newmol(nat,xyz,at,etot,chrg,gen)
     if(.not.allocated(self%mol))then
         self%nmol = 1
         allocate(self%mol(1))
         allocate(self%new(1))
         self%mol(1) = mol
         self%new(1) = .true.
     else
         n = self%nmol + 1
         allocate(dummy(n))
         allocate(btmp(n))
         dummy(1:self%nmol) = self%mol(1:self%nmol)
         dummy(n) = mol
         btmp(1:self%nmol) = self%new(1:self%nmol)
         btmp(n) = .true.
         self%nmol = n
         call move_alloc(btmp, self%new)
         call move_alloc(dummy, self%mol)
     endif
     return
end subroutine ilist_append

  !! QCXMS
 !call fragment_structure(nuc,iat,xyz,3.0_wp,1,0,list)

  ! TAKEN FROM QCxMS
  !  Subroutine for definition of two or more fragments
  !  if at1 = 0 :  global separation in (nonbonded parts), beginning with atom at2
  !  if at1 and at2 <> 0 : define fragments only if a at1-at2 bond (and no other) exists
  !  if at1 and at2 = 0 : delete all fragment assignments
  !  no bond if rcut times the cov.dist.
subroutine fragment_structure(nat,oz,xyz,rcut,at1,at2,frag)
  
     integer  :: at1,at2, nat
     integer  :: i,j
     integer  :: attotal,currentfrag
     integer  :: oz(nat),frag(nat)

     real(wp),intent(in) ::  xyz(3,nat)
     real(wp) :: rcov, r
     real(wp) :: rcut
  
     logical  :: finish
     logical  :: connect(nat,nat)


  

   real(wp), parameter :: autoaa = 0.52917726_wp
   real(wp), parameter :: aatoau = 1.0_wp/autoaa  

   ! Radius used in QCxMS (in au)
   real(wp), parameter :: Rad(118) = aatoau *  [ &
   & 0.32_wp,0.37_wp, & ! H,He
   & 1.30_wp,0.99_wp,0.84_wp,0.75_wp,0.71_wp,0.64_wp,0.60_wp,0.62_wp, & ! Li-Ne
   & 1.60_wp,1.40_wp,1.24_wp,1.14_wp,1.09_wp,1.04_wp,1.00_wp,1.01_wp, & ! Na-Ar
   & 2.00_wp,1.74_wp, & ! K,Ca
   &                 1.59_wp,1.48_wp,1.44_wp,1.30_wp,1.29_wp, & ! Sc-
   &                 1.24_wp,1.18_wp,1.17_wp,1.22_wp,1.20_wp, & ! -Zn
   &                 1.23_wp,1.20_wp,1.20_wp,1.18_wp,1.17_wp,1.16_wp, & ! Ga-Kr
   & 2.15_wp,1.90_wp, & ! Rb,Sr
   &                 1.76_wp,1.64_wp,1.56_wp,1.46_wp,1.38_wp, & ! Y-
   &                 1.36_wp,1.34_wp,1.30_wp,1.36_wp,1.40_wp, & ! -Cd
   &                 1.42_wp,1.40_wp,1.40_wp,1.37_wp,1.36_wp,1.36_wp, & ! In-Xe
   & 2.38_wp,2.06_wp, & ! Cs,Ba
   &         1.94_wp,1.84_wp,1.90_wp,1.88_wp,1.86_wp,1.85_wp,1.83_wp, & ! La-Eu
   &         1.82_wp,1.81_wp,1.80_wp,1.79_wp,1.77_wp,1.77_wp,1.78_wp, & ! Gd-Yb
   &                 1.74_wp,1.64_wp,1.58_wp,1.50_wp,1.41_wp, & ! Lu-
   &                 1.36_wp,1.32_wp,1.30_wp,1.30_wp,1.32_wp, & ! -Hg
   &                 1.44_wp,1.45_wp,1.50_wp,1.42_wp,1.48_wp,1.46_wp, & ! Tl-Rn
   & 2.42_wp,2.11_wp, & ! Fr,Ra
   &         2.01_wp,1.90_wp,1.84_wp,1.83_wp,1.80_wp,1.80_wp,& ! Ac-Pu
   ! from covalent 2009 covalent radii, such that it is complete up to 118
   &                                                         1.49_wp, & ! Am
   &         1.49_wp,1.51_wp,1.51_wp,1.48_wp,1.50_wp,1.56_wp,1.58_wp, & ! Cm-No
   &                 1.45_wp,1.41_wp,1.34_wp,1.29_wp,1.27_wp, & ! Lr-
   &                 1.21_wp,1.16_wp,1.15_wp,1.09_wp,1.22_wp, & ! -Cn
   &                 1.36_wp,1.43_wp,1.46_wp,1.58_wp,1.48_wp,1.57_wp ] ! Nh-Og    
  
     connect(1:nat,1:nat)=.false.
  
     do i = 1,nat-1
        do j = i+1, nat
           r = sqrt((xyz(1,i)-xyz(1,j))**2 + (xyz(2,i)-xyz(2,j))**2 &
           & + (xyz(3,i)-xyz(3,j))**2)
           rcov = rcut * 0.5_wp * (Rad(oz(i)) + Rad(oz(j)))
           if(r.lt.rcov)then
              connect(i,j)=.true.
              connect(j,i)=.true.
           endif
        enddo
     enddo
     if ((at1.eq.0).and.(at2.eq.0)) then
        do i = 1,nat
           frag(i) = 1
        end do
        return
     else
  
        do i = 1,nat
           frag(i) = 0
        end do
  
        frag(at1) = 1
        attotal=1
  
        if (at2.ne.0) then
           connect(at1,at2) = .false.
           connect(at2,at1) = .false.
        endif
  
        finish=.false.
        currentfrag=0
  
        do while (attotal.ne.nat)
  
           currentfrag=currentfrag + 1
  
  !  cycle through atoms and find connected ones
  
           do while (.not.(finish))
              finish=.true.
              do i = 1,nat
                 if (frag(i).eq.currentfrag) then
                    do j = 1,nat
                       if (connect(i,j)) then
                          if (frag(j).eq.0) then
                             frag(j) = currentfrag
                             attotal = attotal + 1
                             finish = .false.
                          elseif (frag(j).eq.currentfrag) then
                             cycle
                          endif
                       endif
                    enddo
                 endif
              enddo
           enddo
  
  ! find the first atom in the next fragment
  
           do i = 1,nat
              if (frag(i).eq.0) then
                 frag(i) = currentfrag + 1
                 attotal = attotal + 1
                 exit
              endif
           enddo
           finish=.false.
        enddo
  
     endif
     return
  end subroutine fragment_structure

end module msmod
