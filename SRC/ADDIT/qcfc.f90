program qcfc

    implicit none
    integer narg 
    integer N,i,j,k,l,Charge,Spin_Mult,rec_char,opt
    integer,allocatable :: MAtom(:)
    real(kind=8),allocatable:: Rn(:,:)
    character(len=80) Ch0,Ch1,Atoms
    character(len=200) fname,sfname
    character(len=2),allocatable:: Atom(:)
    
    narg = IARGC()
    if( narg > 0 ) then
       call getarg(1,fname)
    else
       write(*,*) " input gaussian input file's name. "
       read(*,*) fname
    end if
    
    call Filename(fname,sfname,rec_char)
    open(unit = 10,file = fname)   
    read(10,"(A1)") Ch0
    close(10)

    if ( Ch0 == "#" .or. Ch0 == "%" ) then
        !---------------------------------------------------------------
        open(unit = 20,file = fname)
        
        k = -1
        do while(.true.)
           read(20,"(A1)") Ch1
           k = k + 1
           if( Ch1 == "#" ) exit
        end do
        
        read(20,*)
        read(20,*)
        read(20,*)
        read(20,*)
        
        i = 0
        do while(.true.)
           read(20,"(A2)") Atoms
           i = i + 1
           if( Atoms == "  " ) exit
        end do
        N = i - 1
        allocate(Rn(3,N),Atom(N),MAtom(N))
        close(20) 
        
        open(unit = 30,file = fname)
        do i = 1,k
            read(30,*)
        end do
        
        read(30,*)
        read(30,*)
        read(30,*)
        read(30,*)
        read(30,*) Charge,Spin_Mult
        do i = 1,N
           read(30,*) Atom(i),Rn(1:3,i)
        end do
        close(30)
        
        call Cal_NAtom(N,MAtom,Atom)
        !---------------------------------------------------------------
        write(*,*)
        write(*,*) "Choose you option:"
        write(*,*) "  1: gaussian --> gamess."
        write(*,*) "  2: gaussian --> molpro."
        write(*,*) "  3: gaussian --> nwchem."
        write(*,*) "  4: gaussian --> ORCA."
        read(*,*) opt
        
        if( opt == 1 ) then
            open(unit=100,file = sfname(1:rec_char)//".inp")
            call gamess(N,MAtom,Atom,Rn)
            close(100)
            write(*,*)
            write(*,*) "Done!"
        end if
        
        if( opt == 2 ) then
            open(unit=100,file = sfname(1:rec_char)//".inp")
            call molpro(N,MAtom,Atom,Rn)
            close(100)
            write(*,*)
            write(*,*) "Done!"
        end if
        
        if( opt == 3 ) then
            open(unit=100,file = sfname(1:rec_char)//".nw")
            call nwchem(N,Atom,Rn)
            close(100)
            write(*,*)
            write(*,*) "Done!"
        end if
        
        if( opt == 4 ) then
            open(unit=100,file = sfname(1:rec_char)//".inp")
            call ORCA(N,Atom,Rn)
            close(100)
            write(*,*)
            write(*,*) "Done!"
        end if
                
    else 
       write(*,*) " Wrong input file!"
       stop
    end if
    
end program qcfc
!=======================================================================
subroutine gamess(N,MAtom,Atom,Rn)

    implicit none
    integer i,N,MAtom(N)
    real(kind=8) Rn(3,N)
    character(len=2) Atom(N)
    
    write(100,*) "$BASIS GBASIS=N21 NGAUSS=3 $END"
    write(100,*) "$CONTRL SCFTYP=RHF RUNTYP=ENERGY $END"
    write(100,*)
    write(100,*) "$DATA"
    write(100,"(A5)") "Title"
    write(100,"(A2)") "C1"
    
    do i = 1,N
        write(100,"(A2,F6.1,F20.8,2F14.8)") Atom(i),real(MAtom(i)),Rn(1:3,i)
    end do  
    
    write(100,*) "$END"  
    
end subroutine gamess
!=======================================================================
subroutine molpro(N,MAtom,Atom,Rn)

    implicit none
    integer i,N,MAtom(N)
    real(kind=8) Rn(3,N)
    character(len=2) Atom(N)
    
    write(100,"(A9)") "*** title"
    write(100,*)
    write(100,"(A12)") "gprint,basis"
    write(100,"(A14)") "gprint,orbital"
    write(100,*)
    write(100,"(A12)") "basis, 3-21G"
    write(100,*) 
    write(100,"(A11)") "geomtyp=xyz"
    write(100,"(A10)") "geometry={"
    write(100,"(g0)") N
    write(100,*)
    
    do i = 1,N
        write(100,"(A2,F20.8,2F14.8)") Atom(i),Rn(1:3,i)
    end do  
    
    write(100,"(A1)") "}"  
    write(100,*)
   
    write(100,"(A4)") "{rhf"
    write(100,"(A3,I4,A5)") "wf,",sum(MAtom(:)),",1,0}"  
    write(100,*)
    write(100,"(A3)") "---"
    
end subroutine molpro
!=======================================================================
subroutine nwchem(N,Atom,Rn)

    implicit none
    integer i,N
    real(kind=8) Rn(3,N)
    character(len=2) Atom(N)
    
    write(100,"(A4)") "echo"
    write(100,*)
    write(100,"(A14)") "start molecule"
    write(100,*)
    write(100,"(A13)") 'title "title"'
    write(100,"(A8)") "charge 0"
    write(100,*) 
    write(100,"(A42)") "geometry units angstroms print xyz autosym"
    
    do i = 1,N
        write(100,"(A2,F20.8,2F14.8)") Atom(i),Rn(1:3,i)
    end do  
    
    write(100,"(A3)") "end"
    write(100,*)
    write(100,"(A5)") "basis"    
    write(100,"(A17)") "  * library 3-21G"
    write(100,"(A3)") "end"
    write(100,*)
    write(100,"(A15)") "task scf energy"
    
end subroutine nwchem
!=======================================================================
subroutine ORCA(N,Atom,Rn)

    implicit none
    integer i,N
    real(kind=8) Rn(3,N)
    character(len=2) Atom(N)
    
    write(100,"(A11)") "! RHF 3-21G"
    write(100,"(A8)") "*xyz 0 1"
    
    do i = 1,N
        write(100,"(A2,F20.8,2F14.8)") Atom(i),Rn(1:3,i)
    end do  
    
    write(100,"(A1)") "*"  
    
end subroutine ORCA
!=======================================================================
subroutine Cal_NAtom(M,MAtom,Atom)

    implicit none
    integer M,MAtom(M),i,j
    character(len=2) Atom(M),Atom_Set(54)
    
    Atom_Set = &
    (/'H ',                                                                                'He', &  !  2
      'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &  ! 10
      'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &  ! 18
      'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  ! 36
      'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe'  /) ! 54
    
    do i = 1,M
    
       do j = 1,54
          if( Atom(i) == Atom_Set(j) ) exit
       end do
       MAtom(i) = j
       
    end do
        
end subroutine Cal_NAtom
!=======================================================================
subroutine Filename(fname,sfname,rec_char)

    implicit none
    integer N,i,j,rec_char
    character(len=200) fname,sfname
    
    N = len(fname)
    j = 0
    rec_char = 0
    
    do i = N,1,-1
         if( fname(i:i) == "/" ) then
             j = i
             exit
         end if
    end do 

    do i = j+1,N
         if( fname(i:i) == "." ) exit
         rec_char = rec_char + 1
         sfname(rec_char:rec_char) = fname(i:i)
    end do
       
end subroutine Filename
!=======================================================================
