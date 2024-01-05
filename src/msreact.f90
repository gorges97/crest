!================================================================================!
! This file is part of crest.
!
! Copyright (C) 2020 Stefan Grimme, Philipp Pracht
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

!==============================================================!
! the handler function is the sub-program that is called
! by crest, i.e., within this routine the seperate input file is
! read and the calculation is started.
!==============================================================!
subroutine msreact_handler(env,tim)
    use iso_fortran_env, wp => real64
    use crest_data
    use msmod
    use strucrd
    use zdata ! from JG
    use iomod ! from JG
    use atmasses ! from JG
    implicit none

    type(systemdata) :: env
    type(timer) :: tim

    type(msobj) :: mso  !a new msreact object
    type(coord) :: struc

    interface
        subroutine msreact_topowrap(mol,pair,paths,wboname)
            import :: msmol
            implicit none
            type(msmol) :: mol
            integer :: pair(mol%nat*(mol%nat+1)/2)
            integer :: paths(mol%nat*(mol%nat+1)/2,mol%nat)
            character(len=*),optional :: wboname
        end subroutine msreact_topowrap
    end interface
    integer :: nat
    integer,allocatable :: pair(:)
    integer,allocatable :: paths(:,:)
    integer :: k
    

    !! JC variables
    integer :: nfrags, i, j, ich ! number of current fragmentpairs
    integer :: frag1, frag2, nstruc, incr, ii 
    real(wp) :: mass
    real(wp),allocatable :: xyz(:,:,:) ! xyz coordinates of fragments
    real(wp) :: estart
    integer,allocatable :: at(:) 
    integer ::  nbaseat ! number of lewis basic atoms according to a LMO calculation by xTB
    integer :: fragcount ! number of fragments in structure
    integer, allocatable :: basicatlist(:)
    
    character(len=128) :: atmp
    character(len=128), allocatable :: comments(:)
     
     logical :: equal

    call tim%start(1,'MSREACT')

    !-- read the input coord and put it into the
    !   iso-list as Gen 0 structure
    call struc%open('coord')
    struc%xyz=struc%xyz*bohr !to Angström, from this point on by convention!
    
   
    ! write mass
    mass = molweight(struc%nat,struc%at)
    call wrshort_real("mass",mass)
    call mso%il%append(struc%nat,struc%at,struc%xyz,0.0_wp,env%chrg,0)
    
    
        
    !-- additional input file could be read here
 !   call msinputreader(mso) ! TODO currently via command line arguments

    nat = mso%il%mol(1)%nat
    k=nat*(nat+1)/2
    allocate(pair(k),paths(k,nat))
    call msreact_topowrap(mso%il%mol(1),pair,paths,'wbo')

    ! optimize input structure
    ! necessary?
    call wrxyz('instruc.xyz',struc%nat,struc%at,struc%xyz,'input structure')
    call struc%deallocate()! compare arrays?
    call get_start_energy(env,mso%T,estart) ! TODO fixme remove estart?
  !  write(startcomment,*) "input structure ", estart
  !  call wrxyz('instruc.xyz',struc%nat,struc%at,struc%xyz,trim(startcomment))
 
    write(*,*) "estart is", estart
    !-- do the directory setup and optimizations
    !call msreact(mso,mso%il%mol(1),nat,pair,3)
    
  
        !atrractive potential for H-shifts
      !for H-shifts: xtb -lmo , read lmo.out -> get list of atoms with pi or lone-pair -> H-allowed to shift there
      ! with distance criterium??
      call xtblmo(env,.true.)  
      call readbasicpos(env,nbaseat,basicatlist)

       !-- setting the threads for correct parallelization
      if(env%autothreads)then
        call ompautoset(env%threads,6,env%omp,env%MAXRUN,0) !set global OMP/MKL variable for xtb jobs
      endif
    
    call msreact(env,mso,mso%il%mol(1),nat,pair,nbaseat,basicatlist)
    deallocate(paths,pair)
    !-- read the output and write it to the flist and iso-list
   ! call collect_structures(mso,mso%il,mso%flist)
   
    
    !do some subsequent sorting steps
    call rdensembleparam('MSDIR/products.xyz', nat, nfrags)
    write(*,*)
    write(*,*)"Number of atoms      = ", nat
    write(*,*)"Number of Fragments  = ", nfrags

    ! increase distance between two fragments, this is better for the topology tools(?)
    ! and for the transition state search (we dont want to much interaction between separated fragments)
   ! call fragdistance()
    
    ! first determine center of mass between fragments 
    ! then add to second fragment 3 times the distance vector between both center of masses
     allocate(xyz(3,nat,nfrags),at(nat),comments(nfrags))
    
    call rdensemble('MSDIR/products.xyz',nat,nfrags,at,xyz,comments) 

    ! write optimized starting structure to file
    call move('instruc.xyz','tmpproducts.xyz')
   

    ! increas distance between fragments, better for transistion state search and topology tools? -> makes no difference ....
    call increase_fragdistance(nfrags,xyz,nat,at,mso%fragdist)

      ! wrong etots are written here, changed due to distance !!!!!
      open (newunit=ich,file='MSDIR/products.xyz',status='replace')
    do i = 1,nfrags
      call wrxyz(ich,nat,at,xyz(:,:,i),comments(i))
    end do
    close (ich)
    deallocate(xyz,at,comments)
 
    !include starting structure in cosort, so that if nothing happens with starting structure this one gets sorted out
    ! here, then delete it later from ensemble
    !qcxms2 names every infile fragment.xyz

       !   call copy('MSDIR/products.xyz','products1.xyz')
         ! write(*,*)
        !  write(*,*) "Generated fragment structures written to <products1.xyz>"
    call appendto('MSDIR/products.xyz','tmpproducts.xyz')   

    
     if (env%mstopo) then    
     ! checks additionally if fragment assignment is the same
        call cosort2('tmpproducts.xyz','tmpproducts1.xyz',.false.,.true.)
    elseif (env%msstopo) then
     ! sorts out only according to topovec problem c-c bond rotation changes already topovec
       call sortouttopo('tmpproducts.xyz','tmpproducts1.xyz')
     else 
        call copy('tmpproducts.xyz', 'tmpproducts1.xyz') 
     end if
    
      if(env%msmolbar) then 
       call sortoutmolbar(env,'tmpproducts1.xyz','tmpproducts2.xyz') 
       elseif(env%msinchi) then
      call sortoutinchi('tmpproducts1.xyz','tmpproducts2.xyz')
     else 
     call copy('tmpproducts1.xyz','tmpproducts2.xyz') 
     end if
    ! replace energies by relative energies ?? 
    ! or make ewin dependent on lowest energy? 
    
   
    !call remove('tmpfragments.xyz')
    !call remove('tmpproducts.xyz')
    
    ! sorting step for similar fragments

    ! set energy window for sorting out structures
   ! env%ewin = mso%ewin
   write(*,*) "energy window is", env%ewin
    ! workaround because i dont find where you can rename the input for the cregen routine
  call copy('tmpproducts2.xyz','crest_rotamers_0.xyz') !!! warning if crest_rot
  call remove('crest_rotamers*')
    !!call move('MSDIR/fragments.xyz','crest_rotamers_0.xyz')
    ! in this step structures above en%ewin =10eV and duplicates according to energy or RMSD are sorted out
     !dont do this 
     env%ethr = 0.2_wp ! structures within 0.1 kcal are considered equal 1kcal too much
     !env%rthr = 200_wp ! just set ridicously high to only sort out according to energy
   !  env%trackorigin = .true.   
    call newcregen(env,12)
     
    call copy('crest_conformers.xyz','tmpproducts3.xyz')
    ! comment information lost here ....
    ! only etots are still there


  

    call rdensembleparam('tmpproducts3.xyz', nat, nfrags)
    allocate(xyz(3,nat,nfrags),at(nat),comments(nfrags))
    call rdensemble('tmpproducts3.xyz',nat,nfrags,at,xyz,comments) 
    
   ! sort out input structure from ensemble, detect it by energy
    open (newunit=ich,file='products2.xyz',status='replace')
    do i = 1,nfrags 
      if ( abs(grepenergy(comments(i)) - estart) .lt. 0.0001_wp ) then 
     !if (trim(comments(i))  == "input structure") then
        write(*,*) "remove initial frag", i 
      cycle ! remove initial structure not ideal
      end if
      call wrxyz(ich,nat,at,xyz(:,:,i),comments(i))
    end do
    close (ich)
   write(*,*)
   write(*,*) "Topologically unique structures without initial structure written to file <products2.xyz>" 

   write(*,*)
   write(*,*) "dissociated structures written to <fragments.xyz>"
   write(*,*) "isomers written to <isomers.xyz>"
   
   ! detect number of fragments and write to file
   call detectfragments(env,'products2.xyz')
   
    
   if ( .not. env%mslargeprint) then
    call remove('tmpproducts.xyz')
    call remove('tmpproducts1.xyz')
    call remove('tmpproducts2.xyz')
    call remove('products1.xyz')
    call remove('products2.xyz')
    call remove('crest_conformers.xyz')
    call remove('crest_rotamers_0.xyz')
    call remove('crest_rotamers_1.xyz')
    call remove('crest_best.xyz')
    call remove('gfnff_topo')
    call rmrf('MSDIR')
   end if

  
    ! print only pre-defined number of structures
       nstruc = env%msnfrag
    ! print
     if (nstruc .ne. 0) then
     deallocate(xyz,at,comments)
  call rdensembleparam('products.xyz', nat, nfrags)
    allocate(xyz(3,nat,nfrags),at(nat),comments(nfrags))
    call rdensemble('products.xyz',nat,nfrags,at,xyz,comments) 
   incr = 1
   if (nstruc .ne. 0) then 
   incr = nfrags/nstruc
   if (incr .lt. 1) incr = 1
   write(*,*) "Printing only ",nstruc," selected structures to products.xyz"
    open (newunit=ich,file='products.xyz',status='replace')
    ii = 1
    do i = 1,nfrags, incr    
     
      if (ii .gt. nstruc) exit
      call wrxyz(ich,nat,at,xyz(:,:,i),comments(i))
       ii = ii + 1 
    end do
    close (ich)
    end if
    end if

    
    call write_fragments(env,'products.xyz')
     
    call tim%stop(1)

end subroutine msreact_handler

!==============================================================!
! the main implementation of the msreact algo should go here
!==============================================================!
subroutine msreact(env,mso,mol,nat,pair,nbaseat,basicatlist)
      use iso_fortran_env, wp => real64
      use msmod
      !use crest_data, only : bohr
      use crest_data
      use iomod
      implicit none

      type(msobj) :: mso    !main storage object
      type(msmol) :: mol, molshift  ! xyz etc
      type(systemdata) :: env ! system data

      integer :: nat
      integer :: pair(nat*(nat+1)/2)
      !integer :: paths(nat*(nat+1)/2,nat)
      integer :: lin !this is a function
      integer :: i,j,k, ii
      integer :: p
      integer :: np
      integer :: io

      character(len=:),allocatable :: subdir
      character(len=40) :: pdir
      character(len=512) :: thisdir
      character(len=20) :: gfnver
      real(wp)             :: constr_dist
      real(wp),allocatable :: rcov(:)
      integer, intent(in) :: basicatlist(nbaseat)
      integer :: nbaseat

      allocate(rcov(94))
      call setrcov(rcov)

      !-- main subdirectory handling
      call getcwd(thisdir)
      subdir='MSDIR'
      io = makedir(subdir)
      call chdir(subdir)

    !atrractive potential for H-shifts
    !for H-shifts: xtb -lmo , read lmo.out -> get list of atoms with pi or lone-pair -> H-allowed to shift there
    !with distance criterium??
    if (env%msattrh) then
      write(*,*) "Lewis basic atoms for H-shifts"
      do i = 1, nbaseat
        write(*,*) basicatlist(i)
      end do
    end if
       !-- get specific pairs
      np=0
      do p=1, env%msnbonds    ! bonds in between
!        write(*,'(1x,a,i0,a)') '1,',p+1,' pairs'
         do i=1,nat
           do j=i,nat
              k=lin(i,j)
              if(p.eq.1.and.pair(k).eq.1.and.(mol%at(i).eq.1.or.mol%at(j).eq.1)) cycle
              if(pair(k)==p)then
                 np = np+1 
                 write(pdir,'(a,i0)')'Pair_',np
                 constr_dist = mso%cdist*(rcov(mol%at(i))+rcov(mol%at(j)))*bohr + float(p)
!                write(*,*) mol%at(i),mol%at(j),constr_dist
                 call isodir(mso,trim(pdir),mol,i,j,constr_dist,'rep')
                 ! for  H-shifts attractive potential between H and O or N
               !  if ((mol%at(i) == 1 .and. (mol%at(j) == 7 .or. mol%at(j) == 8)) & 
                ! & .or. (mol%at(i) == 7 .or. mol%at(i) == 8) .and. mol%at(j) == 1) then
                 ! does not really what I want dont understand
                 ! we dont want to get bonded atoms nearer
                 ! only if only one of the atoms is a H-atom
                 ! or take here only O??? .eq. 8 instead of ne 1
                 ! hydrogen has basic atom nearby ?
            !      if (( mol%at(i) .eq. 1 .and. (findloc(basicatlist,j, 1) .ne. 0) ) .or. &
              ! &  ( mol%at(j) .eq. 1 .and. (findloc(basicatlist,i, 1) .ne. 0) )) then

               !  if (p .ge. 2) then
                ! np = np+1 
                 ! write(*,*) "add attractive potential for pair ", np
                ! write(pdir,'(a,i0)')'Pair_',np
                 !constr_dist = 0.25_wp*(rcov(mol%at(i))+rcov(mol%at(j)))*bohr
                 !constr_dist = 0.5_wp*bohr
!                write(*,*) mol%at(i),mol%at(j),constr_dist
                 !call isodir(mso,trim(pdir),mol,i,j,constr_dist)
                ! end if
                !end if
              endif    
           enddo
         enddo
      enddo
    
    if (env%msattrh) then
    write(*,*) "Add attractive potential for H-shifts:"
          !-- attractive potential for H-shifts or distance criterion??
        ! TODO fixme rewrite that only distances are scanned regardless of bonds in between
      do p=1,7   ! bonds in between
!        write(*,'(1x,a,i0,a)') '1,',p+1,' pairs'
         do i=1,nat
           do j=i,nat
              k=lin(i,j)
             if(p.eq.1.and.pair(k).eq.1.and.(mol%at(i).eq.1.or.mol%at(j).eq.1)) cycle
              if(pair(k)==p)then
               
                if (( mol%at(i) .eq. 1 .and. (findloc(basicatlist,j, 1) .ne. 0) ) .or. &
                &  ( mol%at(j) .eq. 1 .and. (findloc(basicatlist,i, 1) .ne. 0) )) then
                  if (p .ge. 2) then
                  ! distance criterion 3.0 Angstroem???
                  ! McLafferty not achieved by this but in two steps made possible with this
                  ! in next fragmentation
                  ! a large distance here accounts for the fact that we often do not 
                 ! have the required conformer present for this rearrangements
                 ! for complicated molecules 5.0 Angstroem necessary, 4.0 at least
                 ! for quinalphos 5 for example, make it as parameter???
                  !  if (sqrt(sum((mol%xyz(:,i)-mol%xyz(:,j))**2)) .lt. 4.0_wp) then
                    if (msmoldist(mol,i,j) .lt. mso%distthr_att) then
                      np = np+1 
                      write(*,*) "add attractive potential for pair ", np," between ",i," and ",j
                      write(pdir,'(a,i0)')'Pair_',np
                      !constr_dist = mso%cdist_att*(rcov(mol%at(i))+rcov(mol%at(j)))*bohr !0.1_wp * 
                      constr_dist = mso%cdist_att*(rcov(mol%at(i))+rcov(mol%at(j)))*bohr + float(p)
                      call isodir(mso,trim(pdir),mol,i,j,constr_dist,'attr')
                    end if
                  end if
                end if
              endif    
           enddo
         enddo
      enddo
    end if
      ! do additional shifting of atoms and subsequent optimization
      write(*,*) "Shift atoms ",env%msnshifts," times and reoptimize"
    
       do ii = 1, env%msnshifts
        call shiftatoms(mso,mol,molshift,ii)
        np = np+1 
        write(pdir,'(a,i0)')'Pair_',np
        call isodiropt(mso,trim(pdir),molshift)
      end do
      ! number of structures explodes here but maybe necessary for planar molecules     

 write(*,*) "Shift atoms",env%msnshifts2," times and apply repulsive potential on bonds of distorted structure"
!nbonds = nbonds - 1
 do ii = 1, env%msnshifts2
  call shiftatoms(mso,mol,molshift,ii)
do p=1,env%msnbonds    ! bonds in between
!        write(*,'(1x,a,i0,a)') '1,',p+1,' pairs'
         do i=1,nat
           do j=i,nat
              k=lin(i,j)
              if(p.eq.1.and.pair(k).eq.1.and.(mol%at(i).eq.1.or.mol%at(j).eq.1)) cycle
              if(pair(k)==p)then
                 np = np+1 
                 write(pdir,'(a,i0)')'Pair_',np
                 constr_dist = mso%cdist*(rcov(molshift%at(i))+rcov(molshift%at(j)))*bohr + float(p)
!                write(*,*) mol%at(i),mol%at(j),constr_dist
                 call isodir(mso,trim(pdir),molshift,i,j,constr_dist,'rep')
                 end if 
           enddo
         enddo
      enddo
enddo

      
      write(*,*) '# of distortions',np
      !-- do the job construction and execution

      call msreact_jobber(np,'Pair_',env%gfnver,.false.)
      call msreact_collect(mol%nat,np,'products.xyz')
      call rename(subdir//'/'//'products.xyz','products.xyz')
      call chdir(thisdir)
      return
end subroutine msreact


!============================================================!
! make a dir for a structure without fragments,
! a controlfile with constraints on atoms A and B (at dist D)
! will be written into the directory
!============================================================!
subroutine isodir(mso,dirname,mol,A,B,D,bias)
    use iso_fortran_env, only : wp => real64
    use msmod
    use iomod
    use strucrd, only : wrxyz
    implicit none
    type(msobj) :: mso
    character(len=*) :: dirname
    character(len=*) :: bias
    type(msmol) :: mol
    integer :: A,B
    real(wp) :: D

    character(len=:),allocatable :: fname
    character(len=20) :: dumm
    integer :: io,ich

    io = makedir(dirname) !create the directory
    
    fname = trim(dirname)//'/'//'struc.xyz'
    open(newunit=ich,file=fname)
    call wrxyz(ich,mol%nat,mol%at,mol%xyz)
    close(ich)

    fname = trim(dirname)//'/'//'.CHRG'
    open(newunit=ich,file=fname)
    write(ich,'(i0)') mol%chrg + 0   ! EI +1, DEA -1, CID 0, give in crest call 
    close(ich)

    fname = trim(dirname)//'/'//'.xc1'
    open(newunit=ich, file=fname)
    write(ich,'(a)') '$scc'
    write(dumm,'(f16.2)') mso%T
    write(ich,'(1x,a,a)')'temp=',adjustl(trim(dumm))
    write(ich,'(a)')'$constrain'
    if (bias == 'rep') write(dumm,'(f16.4)') mso%fc
    if (bias == 'attr') write(dumm,'(f16.4)') mso%fc_attr
    write(ich,'(3x,a,a)')'force constant=',adjustl(trim(dumm))
    write(ich,'(3x,a,1x,i0,a,1x,i0,a,1x,f8.5)') 'distance:',A,',',B,',',D
    close(ich)

    fname = trim(dirname)//'/'//'.xc2'
    open(newunit=ich, file=fname)
    write(ich,'(a)') '$scc'
    write(dumm,'(f16.2)') mso%T
    write(ich,'(1x,a,a)')'temp=',adjustl(trim(dumm))
    write(ich,'(a)') '$opt'
    write(ich,'(1x,a)') 'maxcycle=15'
    write(ich,'(a)') '$write'
    write(ich,'(1x,a)') 'wiberg=true'
    close(ich)

    return
end subroutine isodir

!============================================================!
! make a dir for a structure without fragments,
! a controlfile for an optimization after the atom shifting
! will be written into the directory
!============================================================!
subroutine isodiropt(mso,dirname,mol)
  use iso_fortran_env, only : wp => real64
  use msmod
  use iomod
  use strucrd, only : wrxyz
  implicit none
  type(msobj) :: mso
  character(len=*) :: dirname
  type(msmol) :: mol
  integer :: A,B
  real(wp) :: D

  character(len=:),allocatable :: fname
  character(len=20) :: dumm
  integer :: io,ich

  io = makedir(dirname) !create the directory
  
  fname = trim(dirname)//'/'//'struc.xyz'
  open(newunit=ich,file=fname)
  call wrxyz(ich,mol%nat,mol%at,mol%xyz)
  close(ich)

  fname = trim(dirname)//'/'//'.CHRG'
  open(newunit=ich,file=fname)
  write(ich,'(i0)') mol%chrg + 0   ! EI +1, DEA -1, CID 0, give in crest call 
  close(ich)

  fname = trim(dirname)//'/'//'.xc1'
  open(newunit=ich, file=fname)
  write(ich,'(a)') '$scc'
  write(dumm,'(f16.2)') mso%T
  write(ich,'(1x,a,a)')'temp=',adjustl(trim(dumm))
  close(ich)

  fname = trim(dirname)//'/'//'.xc2' ! will be changed in implementation, we only need to optimize once
  open(newunit=ich, file=fname)
  write(ich,'(a)') '$scc'
  write(dumm,'(f16.2)') mso%T
  write(ich,'(1x,a,a)')'temp=',adjustl(trim(dumm))
  write(ich,'(a)') '$write'
  write(ich,'(1x,a)') 'wiberg=true' ! write wiberg bond orders maybe for later to decect atoms involved in bond changes
  close(ich)
  return
end subroutine isodiropt

!=====================================================================!
! The job construction routine for MSREACT
! (will have to be modified later, for now it is for testing)
!=====================================================================!
subroutine msreact_jobber(ndirs,base,gfnver,niceprint)
     use iso_fortran_env, only : wp => real64
    use msmod
    use iomod
    implicit none
    integer :: ndirs
    character(len=*) :: base
    logical :: niceprint

    character(len=1024) :: jobcall
    character(len=1024) :: jobcall2
    character(len=20) :: gfnver
  integer :: val

    jobcall = ''
    jobcall2 = ''

    write(jobcall,'(a)')  'xtb struc.xyz  --opt loose '//trim(gfnver)//' --input .xc1 > split.out 2>/dev/null'
    write(jobcall2,'(a)') 'xtb xtbopt.xyz --opt crude '//trim(gfnver)//' --input .xc2 > xtb.out 2>/dev/null'
  !  write(jobcall2,'(a)') 'xtb xtbopt.xyz --opt vtight '//trim(gfnver)//' --input .xc2 > xtb.out 2>/dev/null'
    jobcall = trim(jobcall)//' ; '//trim(jobcall2)

    !-- directories must be numbered consecutively
    call opt_OMP_loop(ndirs,base,jobcall,niceprint)
    write(*,*)
    write(*,*) 'done.'
    return
end subroutine msreact_jobber



!=====================================================================!
! A wrapper to generate the topology for a molecule within the
! MSREACT subprogram
!=====================================================================!
subroutine msreact_topowrap(mol,pair,paths,wboname)
    use iso_fortran_env, only : wp => real64
    use msmod
    use zdata
    use crest_data, only : bohr
    implicit none
    type(msmol) :: mol
    integer :: pair(mol%nat*(mol%nat+1)/2)
    !integer :: pair(mol%nat,mol%nat)
    integer :: paths(mol%nat*(mol%nat+1)/2,mol%nat)
    character(len=*),optional :: wboname 
    type(zmolecule) :: zmol

    integer,allocatable :: A(:,:)
    integer,allocatable :: prev(:,:)
    real(wp),allocatable :: E(:,:)
    real(wp),allocatable :: dist(:,:)

    integer :: lpath,i,j,k
    integer :: lin !this is a function
    integer,allocatable :: path(:)
    logical :: ex


    ex=.false.
    if(present(wboname))then
      inquire(file=wboname,exist=ex)
    endif  
    if(ex)then
      call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.false.,wboname)
    else   
       mol%xyz = mol%xyz / bohr !CN based topo requires Bohrs 
       call simpletopo(mol%nat,mol%at,mol%xyz,zmol,.false.,.false.,'')
       mol%xyz = mol%xyz * bohr
    endif 

    allocate(A(mol%nat,mol%nat),E(mol%nat,mol%nat))
    call zadjacent(zmol,A,E)

    allocate(prev(mol%nat,mol%nat),dist(mol%nat,mol%nat))

    call FloydWarshall(mol%nat,A,E,dist,prev)
    allocate(path(mol%nat), source = 0)
    do i=1,mol%nat
      do j=i,mol%nat
       path = 0
       call getPathFW(mol%nat,prev,i,j,path,lpath)
       !write(*,*) path(1:lpath)
       k=lin(i,j)
       pair(k) = lpath - 1 ! number of bonds
       paths(k,:) = path(:)
      enddo
     enddo

    deallocate(dist,prev)
    deallocate(E,A)

    call zmol%deallocate() !clear the zmol memory
    return
end subroutine msreact_topowrap    

!========================================================================!
! collect structures of optimized molecules
! xyz files should still have the same number and order of atoms
!========================================================================!
subroutine msreact_collect(nat,np,outfile)
    use iso_fortran_env, wp => real64
    use strucrd
    use crest_data, only : bohr
    use msmod
    implicit none
    integer :: nat
    integer :: np
    character(len=*) :: outfile
    integer :: ich
    character(len=40) :: pdir
    character(len=:),allocatable :: optfile
    character(len=128) :: newcomment
    integer :: p
    logical :: ex
    integer,allocatable :: at(:)
    real(wp),allocatable :: xyz(:,:)
    real(wp) :: etot


    allocate(at(nat),xyz(3,nat))
    open(newunit=ich,file=outfile)
 
    do p=1,np
      ! write(pdir,'(i0,i0,a,i0)')1,p+1,'Pair_',p2
       write(pdir,'(a,i0)')'Pair_',p
       optfile=trim(pdir)//'/'//'xtbopt.xyz'
       inquire(file=optfile,exist=ex)
       if(ex)then
          call rdcoord(optfile,nat,at,xyz,etot)
          xyz = xyz*bohr
          write(newcomment,'(1x,f18.8,5x,a)')etot,trim(pdir)
          call wrxyz(ich,nat,at,xyz,newcomment)
        end if       
    enddo
    close(ich)
    deallocate(xyz,at)
    return
end subroutine msreact_collect

!============================================================!
! shift atoms of input structure randomly
! to generate more structures
! for planar molecules important 
!============================================================!
subroutine shiftatoms(mso,mol1,mol2,count)
    use iso_fortran_env, only : wp => real64
    use msmod
    use iomod
    use strucrd, only : wrxyz
    implicit none
    type(msobj) :: mso
    type(msmol) :: mol1, mol2
    integer :: natshift
    integer :: natnoh, shiftcount
    real(wp) :: D
    real(wp) :: x, y, shift 
    real(wp) :: dist ! distance two shifts atoms
    character(len=80) :: fname
    character(len=20) :: dumm
    integer :: io,ich, i, j, count, incr
    dist = 0.75_wp ! in angstroem? 0.5 too small, 1.0 maybe too large



    mol2 = mol1
  !  write(*,*) "Number of atoms without H", natnoh
    ! shift random number of atoms randomized selected
  !   natshift = mol1%nat
   !  call Random_Number(x)
    ! natshift = nint(x * mol1%nat)
    ! if (natshift .eq. 0) natshift = 1
   ! write(*,*) "shifted atoms", natshift
   ! incr = natshift/mol1%nat
     
    !if (incr .lt. 1) incr = 1
   ! write(*,*) "increment", incr
    shiftcount = 0
    do i = 1, mol1%nat
     !   if (shiftcount .eq. natshift) exit
         ! do not shift H atoms
       ! if (mol1%at(i) .ne. 1) then        
        do j = 1, 3
            call Random_Number(x)
            ! positive and negative numbers possible
            x = 2.0_wp*x-1.0_wp  
            shift = x* dist 
             mol2%xyz(j,i)=mol1%xyz(j,i) + shift
           
        end do
        shiftcount =  shiftcount + 1
       ! end if
    end do
    ! del
    write(fname,'(a,i0,a)')  'struc',count,'.xyz'
    open(newunit=ich,file=fname)
    call wrxyz(ich,mol2%nat,mol2%at,mol2%xyz)
    close(ich)
    return       
end subroutine shiftatoms

 ! increase distance between fragments, better for transistion state search and topology tools?G:q
subroutine increase_fragdistance(nfrags,xyz,nat,at,fragdist)
    use iso_fortran_env, only : wp => real64
    use msmod
    use iomod
    use strucrd, only : wrxyz
    use axis_module
    implicit none
    integer, intent(in) :: nat
    integer, intent(in) :: at(nat)
    integer, intent(in) :: nfrags
    real(wp), intent(in) :: fragdist
    real(wp), intent(inout) :: xyz(3,nat,nfrags)
    integer, allocatable :: at1(:), at2(:)
    integer, allocatable :: fragi(:)
    integer :: frag1, frag2, i, j, fragcount
    real(wp) :: norm
    real(wp),allocatable :: xyz1(:,:), xyz2(:,:) ! xyz of fragment 1 and 2
    real(wp) :: cmass1(3), cmass2(3) ! center of mass of fragment 1 and 2
    type(msobj) :: mso 
    frag1 = 0
    frag2 = 0
    write(*,'(a,f10.8)') "Increasing distance of fragments by ", fragdist, " angstroem"
    allocate(fragi(nat))
    do i = 1, nfrags 
        call fragment_structure(nat,at,xyz(:,:,i),mso%rcut,1,0,fragi,fragcount) ! better than mrec
          if (fragcount .gt. 1) then 
          do j = 1, nat
            if (fragi(j)==1) then
            frag1 = frag1 +1 
            end if
            if (fragi(j)==2) then
            frag2 = frag2 +1 
            end if
          end do
           allocate(at1(frag1),at2(frag2))
           allocate(xyz1(3,frag1),xyz2(3,frag2))
           frag1 = 0
           frag2 = 0
           do j = 1, nat
            if (fragi(j)==1) then
            frag1 = frag1 +1 
            xyz1(:,frag1)=xyz(:,j,i)
            at1(frag1)=at(j)
            end if
            if (fragi(j)==2) then
            frag2 = frag2 +1 
            xyz2(:,frag2)=xyz(:,j,i)
            at2(frag2)=at(j)
            end if
          end do
          call CMAv(frag1,at1,xyz1,cmass1)
          call CMAv(frag2,at2,xyz2,cmass2)
          norm = sqrt( (cmass2(1) - cmass1(1))**2 + (cmass2(2) - cmass1(2))**2 + (cmass2(3) - cmass1(3))**2)
          do j = 1, nat
            if (fragi(j)==2) then
            xyz(:,j,i) = xyz(:,j,i) + (cmass2-cmass1)/norm * fragdist 
            end if
          end do
          deallocate(at1,at2,xyz1,xyz2)
          end if
    end do

        
end subroutine increase_fragdistance
