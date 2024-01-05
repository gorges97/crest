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
      use strucrd ! from JG
      use iomod
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

          
          real(wp) :: wbormsd = 0.5_wp  !wbo rmsd comp. thr. ??? needed
          real(wp) :: ewin = 500.0_wp ! energy window threshold for fragmentpairs
          real(wp) :: T  = 5000.0_wp ! 5000.0_wp !3000.0_wp !3000.0_wp    electronic temperature! better for fragment generation
          real(wp) :: fc = 0.05_wp      !start fc
           real(wp) :: fc_attr = -0.05_wp      ! fc attractive
          real(wp) :: rcut = 1.3_wp ! cutoff for fragmentation of fragment_structure adapted to 1.3 here
          real(wp) :: cdist = 1.5_wp    !constraing distance scaling factor of rcov
          real(wp) :: cdist_att = 0.5_wp !constraing distance scaling factor of rcov for attractive part
          real(wp) :: distthr_att = 4.0_wp ! distance threshold for attractive part of constraint
          real(wp) :: fragdist = 1.0_wp ! increase distance by 1 angstroem in some cases already too much
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


  ! find fragments first so that we can specify later which fragments should be printed
   subroutine detectfragments(env,fname)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      use atmasses
      use crest_data
      implicit none
      character(len=*) :: fname
      character(len=128),allocatable :: etots(:)
      character(len=512) :: thispath,tmppath1,tmppath2, strucname, tmppath3
      real(wp),allocatable :: xyz(:,:,:)
      real(wp) :: mass
      integer :: nat,nfrags, npoly
      integer :: nc
      integer,allocatable :: at(:), fragi(:)
      integer :: i,r, j, ich, k, ii, ich1, ich2, ich3
      integer :: natf(2)
      integer :: fragcount
      type(systemdata) :: env
      logical :: ex

      npoly = 0
      call rdensembleparam(fname,nat,nfrags)
      allocate(xyz(3,nat,nfrags),at(nat),fragi(nat),etots(nfrags))

      call rdensemble(fname,nat,nfrags,at,xyz,etots) 
      call getcwd(tmppath1)
      ii = 0
    
        open (newunit=ich3,file='products.xyz',status='replace')
      do i=1,nfrags
         call fragment_structure(nat,at,xyz(:,:,i),1.3_wp,1,0,fragi,fragcount) 
         ! do not write structures that were not fragmented, i.e., rearrangements
          if (fragcount .gt. 2) then 
        write(*,*) "More than 2 fragments generated -> Pair ", i," is sorted out " ! we only want to get elementary reactions, so no cascade reactions, but rearrangements are still ok
        npoly = npoly +1
        cycle
         end if

         if(env%msnoiso) then
            if (fragcount .eq. 1) cycle
            ! only write rearranged structures
            elseif(env%msiso) then
               if (fragcount .gt. 1) cycle
         end if  
         ii = ii + 1         
         call wrxyz(ich3,nat,at,xyz(:,:,i),etots(i))

       enddo
       
      if(env%msnoiso) then
         write(*,*) "sorted out ", nfrags - ii,"non-dissociated structures and ", npoly," multiple fragmented"
      elseif(env%msiso) then
         write(*,*) "sorted out ", nfrags - ii,"dissociated structure pairs of which ", npoly," multiple fragmented"
      end if

    end subroutine detectfragments
    subroutine write_fragments(env,fname)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      use atmasses
      use crest_data
      implicit none
      character(len=*) :: fname
      character(len=128),allocatable :: etots(:)
      character(len=512) :: thispath,tmppath1,tmppath2, strucname, tmppath3
      character(len=40) :: sumform, sumformula
      real(wp),allocatable :: xyz(:,:,:)
      real(wp) :: mass
      integer :: nat,nfrags, npoly
      integer :: nisomers, nfragpairs
      integer :: nc
      integer,allocatable :: at(:), fragi(:)
      integer :: i,r, j, ich, k, npairs, ich1, ich2, ich3, ii
      integer , allocatable :: atf(:,:) ! atomtypes of fragment
      integer :: natf(2)
      integer :: fragcount ! number of fragments of structure
      type(systemdata) :: env
      type(msobj) :: mso 
      logical :: ex


      nisomers = 0
      nfragpairs = 0
      npoly = 0
      call rdensembleparam(fname,nat,nfrags)
      allocate(xyz(3,nat,nfrags),at(nat),fragi(nat),etots(nfrags))
      allocate(atf(2,nat)) 
      call rdensemble(fname,nat,nfrags,at,xyz,etots) 
      call getcwd(tmppath1)
      npairs = 0
       write(*,*) "writing fragments to <products.xyz> and directories:"
      
       write(*,'(a)')'========================================================'
       write(*,'(a,i0)') " number of generated structures: ", nfrags
       write(*,*) " directory | fragment type | energy"
      open (newunit=ich1,file='isomers.xyz',status='replace')
       open (newunit=ich2,file='fragmentpairs.xyz',status='replace')
        open (newunit=ich3,file='products.xyz',status='replace')
      do i=1,nfrags
          call fragment_structure(nat,at,xyz(:,:,i),mso%rcut,1,0,fragi,fragcount) 
         ! do not write structures that were not fragmented, i.e., rearrangements
          if (count(fragi==3) .gt. 0) then 
        ! write(*,*) "More than 2 fragments generated -> Pair ", i," is sorted out "
         npoly = npoly +1
        cycle
         end if
         !first write isomers
         if (fragcount .eq. 1) then 
            call wrxyz(ich1,nat,at,xyz(:,:,i),etots(i))
         end if
         ! and write fragments
         if (fragcount .gt. 1) then 
            call wrxyz(ich2,nat,at,xyz(:,:,i),etots(i))
         end if

         if(env%msnoiso) then
            if (fragcount .eq. 1) then 
                cycle
            end if 
            ! only write rearranged structures
        elseif(env%msiso) then
            if (fragcount .gt. 1) then 
                cycle
            end if 
         end if  
         npairs = npairs + 1         
         write(tmppath2,'(a,i0)')"p",npairs
         r = makedir(trim(tmppath2))
         ! write to output
         ! check if fragment or isomer
          ii=count(fragi==2)
            if (ii.ne.0) then
               write(*,'(a4,a20,a20)') trim(tmppath2),"fragmentpair",trim(etots(i))
               strucname='pair.xyz'
                nfragpairs = nfragpairs + 1
            else 
               write(*,'(a4,a20,a20)') trim(tmppath2),"isomer",trim(etots(i))
               strucname='isomer.xyz'
               nisomers = nisomers + 1
            end if
        
         call chdir(tmppath2)
         call wrxyz(trim(strucname),nat,at,xyz(:,:,i),etots(i))
         ! just write mass for pairs too
         mass = 0
         do k = 1 , nat
            mass = mass + ams(at(k))
         end do 
         call wrshort_real("mass",mass)  
         call wrxyz(ich3,nat,at,xyz(:,:,i),etots(i))
         call chdir("..")
         if (.not. env%msiso .and. ii .ne. 0) then
         do j = 1, 2 ! currently hardcoded only 2 fragments allowed 
            
            natf(j) = count(fragi==j)
            if (natf(j).ne.0) then
            write(tmppath3,'(a,i0)') trim(tmppath2)//"f",j
            r = makedir(trim(tmppath3))
            call chdir(tmppath3)
            !
            !write(strucname,'(a,i1,a,i1,a)') 'fragment', i,'-', j,'.xyz'
            strucname='fragment.xyz'
            open (newunit=ich,file=strucname,status='replace')
            write(ich,*) natf(j)
            write(ich,*)
            mass = 0
            atf(j,:) = 0
            do k = 1,nat
                if(fragi(k) == j) then !only write xyz if fragment really exists
                  write(ich,'(a2,5x,3F18.8)')  i2e(at(k)),xyz(1:3,k,i) 
                  mass = mass + ams(at(k))
                  atf(j,k) = at(k)
                end if
            end do
            close (ich)
             sumformula = sumform(nat,atf(j,:))
                  ! write also atomic mass here
             call wrshort_real("mass",mass)      
           !   write(*,*) " directory | sumformula | mass"
              write(*,'(a6,i0,5x,a20,9x,f9.5)') trim(tmppath2)//"f",j,trim(sumformula),mass
             ! write info if fragment is isomer
         !    if (count(fragi==2) .eq. 0)   call touch("isomer")  
            end if
            call chdir("..")
         end do
         end if
         call chdir(tmppath1)
       enddo
      call chdir(thispath)
      call wrshort_int('npairs',npairs)
      if(env%msnoiso) then
     ! write(*,*) "sorted out ", nfrags - npairs,"non-dissociated structures and ", npoly," multiple fragmented"
      elseif(env%msiso) then
     ! write(*,*) "sorted out ", nfrags - npairs,"dissociated structure pairs of which ", npoly," multiple fragmented"
      end if
      write(*,'(a)')'========================================================'
     write(*,*)
     write(*,*) "Number of generated isomers: ", nisomers
     write(*,*) "Number of generated fragmentpairs: ", nfragpairs
    end subroutine write_fragments
    ! sort out duplicates according to inchi with obabel
    ! not really necessary because even for LSD only 5 duplicates removed by this
   subroutine sortoutinchi(fin,fout)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      use atmasses
      use crest_data
      implicit none
      type(msobj) :: mso 
      integer :: nat, nfrags
      integer :: io, ich, i, j
      character(len=*) :: fin
      character(len=*) :: fout
      character(len=128),allocatable :: comment(:)
      character(len=512) :: jobcall
      character(len=1024), allocatable :: inchis(:)
      integer,allocatable :: at(:), fragi1(:), fragi2(:)
      integer :: fragcount1, fragcount2
      real(wp),allocatable :: xyz(:,:,:)
      logical, allocatable :: double(:)
      logical :: equal
      write(*,*)
      write(*,*) "Calling obabel for sorting out duplicates by inchicode"
      write(jobcall,'(a)') 'obabel -i xyz '//trim(fin)//' -o inchi > inchicodes 2>obabelout'
      call rdensembleparam(fin,nat,nfrags)
      allocate(xyz(3,nat,nfrags),at(nat),comment(nfrags),fragi1(nfrags),fragi2(nfrags))
      call rdensemble(fin,nat,nfrags,at,xyz,comment)
      
      call execute_command_line(trim(jobcall))
      open (newunit=ich,file='inchicodes',status="old", action="read")
      allocate(inchis(nfrags),double(nfrags))
      
      do i = 1, nfrags
         read(ich,"(a)") inchis(i)
      end do
      close(ich)
      double = .false.
      do i = 1, nfrags
         if(double(i)) cycle
       !  call fragment_structure(nat,at,xyz(:,:,i),mso%rcut,1,0,fragi1,fragcount1)
         do j = i+1, nfrags
            if (.not. double(j)) then 
               ! call fragment_structure(nat,at,xyz(:,:,i),mso%rcut,1,0,fragi2,fragcount2)
           ! call arrcomp(nat,fragi1,nat,fragi2,equal) ! sometimes inchis falsely detects non-bonded structures as bonded -> can throw out too much  
           !if(.not. equal) cycle ! too strict criterion only equal for same atom ordering
            !at least count if numbers add up for fragments
            ! well this is bullshit 
              !  if (count(fragi1==1) .eq. count(fragi2==1) .and. count(fragi1==2) .eq. count(fragi2==2 ) & 
               ! & .and. inchis(i) == inchis(j) ) then  
                 if (inchis(i) == inchis(j)) then   
                    write(*,*) j," is a duplicate of ", i 
                    double(j) = .true.
                else 
                    double(j) = .false.
                end if
            end if
         end do
      end do

      open (newunit=ich,file=fout)
      j = 0
    do i = 1,nfrags
      if ( .not. double(i)) then 
        call wrxyz(ich,nat,at,xyz(:,:,i),comment(i))
      else 
     
        j = j + 1
      end if
    end do
    close (ich)
   write(*,*) "sorted out ", j," fragments accoding to inchi"
   write(*,*) "number of fragmens is now:", nfrags - j," written to ", trim(fout)

   end subroutine sortoutinchi

   subroutine sortouttopo(fin,fout)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      use atmasses
      use crest_data
      use zdata
      implicit none
      integer :: nat, nfrags
      integer :: io, ich, i, j
      character(len=*) :: fin
      character(len=*) :: fout
      character(len=128),allocatable :: comment(:)
      character(len=512) :: jobcall
      integer,allocatable :: at(:), topo1(:), topo2(:)
      integer :: ntopo
      real(wp),allocatable :: xyz(:,:,:)
      logical, allocatable :: double(:)
      logical :: equal


      call rdensembleparam(fin,nat,nfrags)
      allocate(xyz(3,nat,nfrags),at(nat),comment(nfrags))
      call rdensemble(fin,nat,nfrags,at,xyz,comment) 

      allocate(double(nfrags))
      ntopo = nat*(nat+1)/2
      allocate(topo1(ntopo),topo2(ntopo))
      
   
      double = .false.
      do i = 1, nfrags
         if(double(i)) cycle
         call quicktopo(nat,at,xyz(:,:,i),ntopo,topo1)
         do j = i+1, nfrags
            if (.not. double(j)) then
                call quicktopo(nat,at,xyz(:,:,j),ntopo,topo2)
                call arrcomp(ntopo,topo1,ntopo,topo2,equal)
           write(*,*) i 
            write(*,*) topo1
            write(*,*) j 
           write(*,*) topo2
                if (equal) then
                    write(*,*) j," is a duplicate of ", i 
                    double(j) = .true.
                else 
                    double(j) = .false.
                end if
            end if
         end do
      end do

      open (newunit=ich,file=fout)
      j = 0
    do i = 1,nfrags
      if ( .not. double(i)) then 
      call wrxyz(ich,nat,at,xyz(:,:,i),comment(i))
      else 
     
        j = j + 1
      end if
    end do
    close (ich)
   write(*,*) "sorted out ", j," fragments accoding to topology"
   write(*,*) "number of fragmens is now:", nfrags - j," written to ", trim(fout)

   end subroutine sortouttopo

    ! sort out duplicates according to molbar
    ! molbar standalone needed in path
   subroutine sortoutmolbar(env,fin,fout)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      use atmasses
      use crest_data
      implicit none
      integer :: nat, nfrags
      integer :: io, ich, i, j, iocheck
      integer :: nunique
      character(len=*) :: fin
      character(len=*) :: fout
      character(len=80) :: fname
      character(len=512), allocatable :: barcodes(:)
      character(len=128), allocatable :: comments(:)

      character(len=512) :: jobcall
      integer,allocatable :: at(:)
      real(wp),allocatable :: xyz(:,:,:)
      logical, allocatable :: double(:)
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA

      write(*,*)
      write(*,*) "Calling molbar for sorting out duplicates by molbar"
     
      call rdensembleparam(fin,nat,nfrags)
      allocate(xyz(3,nat,nfrags),at(nat),comments(nfrags))
      call rdensemble(fin,nat,nfrags,at,xyz,comments) 
      write(*,*) "calling molbar with ",env%threads," cores"
      do i = 1, nfrags 
        write(fname,'(i0,a)') i,"temp.xyz"
        open (newunit=ich,file=fname)  
        call wrxyz(ich,nat,at,xyz(:,:,i),comments(i))
        close(ich)    
      end do
      ! TAKE CARE, molbar does not print barcodes in the same order as the structures are given
    !  write(jobcall,'(a,i0,a)') 'molbar *temp.xyz -s -T ', env%threads, ' > molbar.out'
    ! correct conda environment can be hard to set up on clusters so that parallelization of molbar works
      write(jobcall,'(a)') 'molbar *temp.xyz -s  > molbar.out'
      call execute_command_line(trim(jobcall))
      
      !identify and remove duplicates
      allocate(double(nfrags),barcodes(nfrags))
     
      ! read lines
      do i = 1, nfrags
        write(fname,'(i0,a)') i,'temp.mb'
        open (newunit=ich,file=fname,status="old", action="read")
        read(ich,'(a)',iostat=iocheck) barcodes(i)
        close(ich)
      end do
      
      ! identify duplicates
      double = .false.
    do i = 1, nfrags
        if (double(i)) cycle
        do j = i + 1, nfrags
            if (.not. double(j)) then
                if (trim(barcodes(i)) == trim(barcodes(j))) then
                    double(j) = .true.
                    write(*,*) j," is a duplicate of ", i
                  !  write(*,*) trim(barcodes(i))
                   ! write(*,*) trim(barcodes(j))
                else
                cycle
                end if
            end if
        end do
    enddo
    close(ich)
    nunique = 0
    open (newunit=ich,file=fout)
    do i = 1, nfrags
        if (.not. double(i)) then
            call wrxyz(ich,nat,at,xyz(:,:,i),comments(i)) 
            nunique = nunique  + 1
        end if
    end do
    close (ich)
   write(*,*) "sorted out ", nfrags - nunique," fragments accoding to molbar"
   write(*,*) "number of fragmens is now:", nunique," written to ", trim(fout)
   !remove some files
   do i = 1, nfrags
      write(fname,'(i0,a)') i,'temp.xyz'
      call remove(trim(fname))
      write(fname,'(i0,a)') i,'temp.mb'
      call remove(trim(fname))
   end do

   end subroutine sortoutmolbar

! DEPRECATED AS NEW MOLBAR DOESNT HAVE THIS FUNCTIONALITY ANYMORE
   ! sort out duplicates according to molbar
    ! molbar standalone needed in path
   subroutine oldsortoutmolbar(env,fin,fout)
      use iso_fortran_env, wp => real64
      use iomod
      use strucrd, only: rdensembleparam,rdensemble,wrxyz
      use atmasses
      use crest_data
      implicit none
      integer :: nat, nfrags
      integer :: io, ich, i, j, iocheck, k
      integer :: nunique, nn
      character(len=*) :: fin
      character(len=*) :: fout
      character(len=80) :: ftemp, line, fname
      real(wp),allocatable :: en(:)
      character(len=512) :: jobcall
      character(len=1024), allocatable :: inchis(:)
      integer,allocatable :: at(:), ind(:)
      real(wp),allocatable :: xyz(:,:,:)
      real(wp) :: xx(100)
      logical, allocatable :: double(:)
      type(systemdata) :: env    ! MAIN STORAGE OS SYSTEM DATA

      write(*,*)
      write(*,*) "Calling molbar for sorting out duplicates by inchicode"
     
      call rdensembleparam(fin,nat,nfrags)
      allocate(xyz(3,nat,nfrags),at(nat),en(nfrags))
      call rdensemble(fin,nat,nfrags,at,xyz,en) 
      write(*,*) "calling molbar with ",env%threads," cores"
      do i = 1, nfrags 
         write(ftemp,'(i0,a)') i,"temp.xyz"
          open (newunit=ich,file=ftemp)  
          call wrxyz(ich,nat,at,xyz(:,:,i),en(i))
          close(ich)
          
      end do
      write(jobcall,'(a,i0)') 'molbar *temp.xyz -s -T ', env%threads
          call execute_command_line(trim(jobcall))
      !remove duplicates with molbar files written to unique.csv
      write(jobcall,'(a)') 'molbar *.mb -r'
      call execute_command_line(trim(jobcall))

      open (newunit=ich,file='unique.csv',status="old", action="read")
      ! first line just en
      nunique = -1
      do
         read(ich,'(a)',iostat=iocheck)
         if (iocheck < 0) exit 
            nunique = nunique + 1
      enddo
      rewind(ich)
      allocate(ind(nunique))
      read(ich,*) line
      do i = 1, nunique
         read(ich,'(a)',iostat=iocheck)line
         if (iocheck < 0) exit 
         do k=1,80
            if(line(k:k) == ',')line(k:k)=' '
          enddo
         call readl(line,xx,nn)
         ind(i) = xx(2)
      enddo
      close(ich)

      open (newunit=ich,file=fout)
    do i = 1,nunique
      call wrxyz(ich,nat,at,xyz(:,:,ind(i)),en(ind(i))) 
    end do
    close (ich)
   write(*,*) "sorted out ", nfrags - nunique," fragments accoding to molbar"
   write(*,*) "number of fragmens is now:", nunique," written to ", trim(fout)
   !remove some files
   do i = 1, nfrags
      write(fname,'(i0,a)') i,'temp.xyz'
      call remove(trim(fname))
      write(fname,'(i0,a)') i,'temp.mb'
      call remove(trim(fname))
   end do
   end subroutine oldsortoutmolbar


   ! TODO Modify this to read more advanced input
   ! maybe with toml
   subroutine msinputreader(mso)
    type(msobj) :: mso 

    logical :: ex
     inquire(file="sumreac",exist=ex)
     if (ex) then
        call rdshort_real('sumreac',mso%ewin)
     end if
   end subroutine msinputreader


   !! --------------------------------------------------------------------------------------
!  Sort out all topologically equivalent structures (i.e. conformers) adjusted for NCI structures
!! !!TODO not ideal, because for example for homotopic hydrogens it now matters again
! which is dissociated
!--------------------------------------------------------------------------------------
subroutine cosort2(iname,oname,wrscoord,verbose)
      use iso_fortran_env, only : wp => real64
      use iomod
      use strucrd, only: wrc0,rdensembleparam,rdensemble,wrxyz
      use crest_data, only: bohr
      implicit none

      type(msobj) :: mso   
      character(len=*),intent(in) :: iname
      character(len=*),intent(in) :: oname
      logical,intent(in)          :: wrscoord

      integer :: j,k,q,p,r
      real(wp),allocatable :: xyz(:,:,:),xyztmp(:,:)
      character(len=128), allocatable :: comments(:)
      real(wp),allocatable :: cn(:),bond(:,:)
      integer,allocatable :: at(:),group(:), fragi1(:), fragi2(:)
      real(wp) :: dE
      character(len=10),allocatable :: itensr(:),itensl(:)   !identifier tensor
      character(len=80) :: str
      integer :: n,nall,nonh,gc,sgc,tgc
      integer :: fragcount1, fragcount2

      logical :: verbose
      logical :: equal
      integer :: ochan,ich

      write(*,*)
      write(*,*)'==================================================='
      write(*,'(a)')' Identifying topologically equivalent structures:'

      call rdensembleparam(iname,n,nall)
      allocate(xyz(3,n,nall),comments(nall),group(nall),at(n),xyztmp(3,n),cn(n),bond(n,n),fragi1(n),fragi2(n))
      call rdensemble(iname,n,nall,at,xyz,comments)
     
      call countnonh(n,at,nonh)
      !allocate(itens(nonh,nall))
      allocate(itensr(nonh),itensl(nonh))

!---- identifier sorting loops
      gc=1
      group=0
   
      do r=1,nall
        if(group(r).ne.0)cycle
        xyztmp(:,:)=xyz(:,:,r)/bohr
       ! call fragment_structure(n,at,xyz(:,:,r),mso%rcut,1,0,fragi1,fragcount1)
        call get_itens(n,xyztmp,at,nonh,itensr)
        do q=1,nall
          if(r.eq.q)then
            group(r)=gc
            tgc=1
            cycle
          endif
          xyztmp(:,:)=xyz(:,:,q)/bohr
          !call fragment_structure(n,at,xyz(:,:,q),mso%rcut,1,0,fragi2,fragcount2)
          !call arrcomp(n,fragi1,n,fragi2,equal)
          ! if(.not. equal) cycle ! too strict criterion
          ! at least count if numbers add up for fragments
          ! bullshit
        !   if (count(fragi1==1) .ne. count(fragi2==1) .or. count(fragi1==2) .ne. count(fragi2==2 )) cycle  
         ! different fragmented
          call get_itens(n,xyztmp,at,nonh,itensl)
          sgc=0
          write(*,*) "compare ",r," with ",q
          do p=1,nonh
             !if(itens(p,r)==itens(p,q))then
             !!if (fragi1(p)==fragi2(p)) then
             if(itensr(p)==itensl(p))then
             write(*,*) itensr(p),itensl(p)
               sgc=sgc+1
             else
               exit
            endif
            !!else 
            !!   exit
            !!endif
          enddo
          if(sgc.eq.nonh)then
             group(q)=gc
             tgc=tgc+1
             write(*,*) "structure ", q," is equivalent to ", r
          endif
        enddo
        if(tgc.gt.1)then
        write(*,'(a,i0,a,i0,a)')' Equivalent to ',r,'. structure: ' &
        &                        ,tgc,' structure(s).'
        endif
        gc=gc+1        
      enddo
      write(*,'(a)')' Done.'
      write(*,'(a,a,a)')' Appending file <',trim(oname),'> with structures.'
      write(*,*)

      open(newunit=ochan,file='tmp')

      write(ochan,*)'==================================================='
      write(ochan,*)'============= ordered structure list =============='
      write(ochan,*)'==================================================='
      if(wrscoord)then
      write(ochan,'(a,a)')' written to <scoord.*> and ',trim(oname)
      write(ochan,*)
      write(ochan,'('' scoord.*     ΔE(kcal/mol)   Etot(Eh)  origin'')')
      else
      write(ochan,'(a,a,a)')' written to file <',trim(oname),'>'
      write(ochan,*)
      write(ochan,'('' structure    ΔE(kcal/mol)   Etot(Eh)  origin'')')
      endif

      open(newunit=ich,file=oname)
      k=0
      do j=1,gc-1
         do r=1,nall
            if(group(r).eq.j)then
              k=k+1
              dE=(grepenergy(comments(r))-grepenergy(comments(1)))*627.5095_wp
              write(ochan,'(i5,6x,F10.2,4x,a)') &
              & k,dE,comments(r)
              call wrxyz(ich,n,at,xyz(:,:,r),comments(r))
              where(group.eq.j)group=0
              if(wrscoord)then  ! write new scoord.* if necessary
                write(str,'(''scoord.'',i0)')k
                call wrc0(trim(str),n,at,xyz(:,:,r))
              endif
            else
              cycle
            endif
         enddo
      enddo
      close(ich)
      close(ochan)

      gc=gc-1

      if(nall.ne.gc)then
        write(*,'(a,i0,a,a,a)')' Initial ',nall,' structures from file ',trim(iname),' have'
        write(*,'(a,i0,a)')' been reduced to ',gc,' topologically unique structures.'
      else
        write(*,'(a,i0,a,a,a)')' All initial ',nall,' structures from file ',trim(iname),' are unique.'
      endif

      !deallocate(itens)
      deallocate(itensl,itensr)
      deallocate(bond,cn,xyztmp,at,comments,xyz)
      
      if(verbose) call cat('tmp')
      call remove('tmp')
      !call remove('identify')

end subroutine cosort2

!--------------------------------------------------------------------------------------------
! A quick xtb geometry optimization in xyz coordinates to get starting energy
!--------------------------------------------------------------------------------------------
subroutine get_start_energy(env,etemp,etot)
         use iso_fortran_env, only : wp => real64
         use iomod
         use crest_data
         use strucrd
         implicit none
         type(systemdata) :: env
         character(len=80) :: fname,pipe
         character(len=512) :: jobcall
         logical :: fin
         character(len=256) :: atmp
         integer :: ich,iost,io,i
         type(coord) :: mol
         integer :: ntopo
         integer,allocatable :: topo(:)
         real(wp), intent(out) :: etot
         real(wp) :: etemp ! electronic temperature in K
         logical :: tchange = .false.
         logical :: ldum

!---- small header
         write(*,*)
         call smallhead('xTB Geometry Optimization')
!---- some options
         pipe=' > xtb.out 2>/dev/null'
         call remove('gfnff_topo')
         if(.not.env%chargesfile)call remove('charges')
         call remove('grad')
         call remove('mos')
         call remove('xtbopt.log')
         call remove('xtbrestart')

          fname = '.CHRG'
          open(newunit=ich,file=fname)
          write(ich,'(i0)') env%chrg  ! EI +1, DEA -1, CID 0, give in crest call 
          close(ich)

!---- input xyz file
         fname='instruc.xyz'
      
         write(jobcall,'(a,1x,a,f10.4,1x,a,1x,a,1x,a)') &
         &     trim(env%ProgName),trim(fname)//" --opt vtight --etemp ",etemp,trim(env%gfnver),trim(env%solv),trim(pipe)
         
        
         call execute_command_line(trim(jobcall), exitstat=io)

         call minigrep('xtb.out','optimized geometry written to:',fin)
         if(.not.fin)then
             write(*,*)
             write(*,*) ' Initial geometry optimization failed!'
             write(*,*) ' Please check your input.'
             error stop
         endif
         write(*,*) 'Geometry successfully optimized.'
!----rename optimized file back     
         call move('xtbopt.xyz','instruc.xyz')
        
         call grepval('xtb.out',"| TOTAL ENERGY",ldum,etot)
!---- cleanup
         call remove('xtb.out')
         call remove('energy')
         if(.not.env%chargesfile)call remove('charges')
         call remove('grad')
         call remove('mos')
         call remove('xtbopt.log')
         call remove('xtbtopo.mol')
          call remove('.xtbopttok')
         call remove('xtbrestart')
         call remove('gfnff_topo')
end subroutine get_start_energy


! this subroutine reads the lmo.out file of an xtb -lmo calculation
! and identifies the atoms which have a LP or pi-Orbital
subroutine readbasicpos(env, nbaseat, basicatlist)
        use iso_fortran_env, only : wp => real64
        use iomod
        use crest_data
        implicit none
        type(systemdata) :: env
        integer :: i, io, ich, at1, at2, j
        integer, intent(out) ::  nbaseat
        integer, intent(out), allocatable :: basicatlist(:)
        integer, allocatable :: dumlist(:)
        character(len=512) :: tmp
        character(len=64) :: fname
        character(len=64) :: type, dumc
        real(wp) :: dumr
        logical :: ex
        fname='lmo.out'
        inquire(file=fname,exist=ex)
        if (.not. ex) then 
           write(*,*) "lmo.out not found"
           stop
        end if
        nbaseat = 0
        open(newunit=ich,file=fname)
        do
           read(ich,'(a)',iostat=io) tmp
           if (index(tmp,'files:') .ne. 0) exit
           if (index(tmp,'pi ') .ne. 0)  nbaseat = nbaseat + 2 ! count first two  atoms of pi or delpi bond with highest participation and ignore the rest
           if (index(tmp,'LP ') .ne. 0)  nbaseat = nbaseat + 1 ! count LP
        end do
        allocate(dumlist(nbaseat))
        dumlist = 0
        j = 0
        rewind(ich)
        do
           read(ich,'(a)',iostat=io) tmp
           if (index(tmp,'starting deloc pi') .ne. 0) exit
           if (index(tmp,'files:') .ne. 0) exit
           if (index(tmp,'pi ') .ne. 0)  then 
            write(*,*) trim(tmp)
               backspace(ich)
               if (tmp(64:64) == ' ') then ! if element symbol has one character, there is a space and we use this routine
                   if (tmp(80:80) == ' ') then
                    read(ich,*) dumr, type, dumr, dumr,  dumr, dumr, dumr, at1, dumc, dumc , dumr, at2, dumc, dumc, dumr ! first element has one and second has one character
                   else
                     read(ich,*) dumr, type, dumr, dumr,  dumr, dumr, dumr, at1, dumc, dumc , dumr, at2, dumc, dumr ! first element has one and second has two character
                   end if 
               else
                   if (tmp(80:80) == ' ') then 
                       read(ich,*) dumr, type, dumr, dumr,  dumr, dumr, dumr, at1, dumc, dumr, at2, dumc, dumc, dumr  ! first element has two and second has pm character
                   else 
                       read(ich,*) dumr, type, dumr, dumr,  dumr, dumr, dumr, at1, dumc, dumr, at2, dumc, dumr ! first element has two and second has two character
                   end if
               end if
               write(*,*) "at1 and at2 are", at1, at2
               if(findloc(dumlist, at1,1) .eq. 0) then 
                   j= j + 1
                   dumlist(j) = at1
               end if
               if(findloc(dumlist, at2,1) .eq. 0) then 
                   j= j + 1
                   dumlist(j) = at2
               end if
           end if

           if (index(tmp,'LP ') .ne. 0)  then 
               backspace(ich)
               if (tmp(64:64) == ' ') then ! if element symbol has one character, there is a space and we use this routine
                   read(ich,*) dumr, type, dumr, dumr,  dumr, dumr, dumr, at1, dumc, dumc , dumr
               else ! if element symbol has two characters, there is no space and we use this routine
                   read(ich,*) dumr, type, dumr, dumr,  dumr, dumr, dumr, at1, dumc , dumr
               end if
               if(findloc(dumlist, at1,1) .eq. 0) then 
                   j= j + 1
                   dumlist(j) = at1
               end if
           end if
        end do
        close(ich)
        basicatlist = pack(dumlist, dumlist .ne. 0) ! sort out zeroes
        nbaseat=size(basicatlist)
        ! call remove('lmo.out')
end subroutine readbasicpos

  !! QCXMS
 !call fragment_structure(nuc,iat,xyz,3.0_wp,1,0,list)

  ! TAKEN FROM QCxMS
  !  Subroutine for definition of two or more fragments
  !  if at1 = 0 :  global separation in (nonbonded parts), beginning with atom at2
  !  if at1 and at2 <> 0 : define fragments only if a at1-at2 bond (and no other) exists
  !  if at1 and at2 = 0 : delete all fragment assignments
  !  no bond if rcut times the cov.dist.
subroutine fragment_structure(nat,oz,xyz,rcut,at1,at2,frag,fragcount) ! better than mrec
  
    integer  :: at1,at2, nat
    integer  :: i,j
    integer  :: attotal,currentfrag
    integer  :: oz(nat),frag(nat)
    integer  :: fragcount  
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

    ! cycle through atoms and find connected ones

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
    
    do i = 1, 3 ! is enough, we only need to know we have 1, 2 or more than 2 fragments
        if (count(frag==i) .gt. 0) then 
            fragcount = i
        end if
    enddo
    return
end subroutine fragment_structure
end module msmod
