program main
    
    call Amesp
   
end program main
!=======================================================================
subroutine Amesp    
     !Amesp : Amateurish molecular electronic structure program.(for linux)
     !Copyright (C) 2016 ; Author: (YingFeng Zhang)Warm_Cloud <3355196386@qq.com>
     !Begin from 2016-03-09. 
!-----------------------------------------------------------------------
!
!    N : Size of the matrix ; M : Number of atoms.
!    Nmem : If N < Nmem ,store Two_Electron in memory.
!    Neri : Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2),Neri=N(N<=Nmem),Neri=1(N>Nmem)
!    Ee : Energy of electrons ; E_Nuc : Energy of nucleus.
!    Me : Nucleus Charge.
!    MAtom : Atom index ; NAtom : Atom index couple with basis sets.
!    rec_char : Length of file name with out extension.
!    Rn : Coordinate of atmos . Rt : Coordinate of atmos couple with basis sets.
!    bas_num : Dimension of every basis set.
!    MoCu(:,i) : Wave function . En : Energy of orbit .
!    c_Basis,a_Basis,l_Basis : c_Basis*x^l_Basis*Exp[-a_Basis*r^2].
!    H_Core : One-electron integral ; Two_Electron : Two-electron integral .
!    f_Mo : <i|T + 1/R|j> ; g_Mo : <i;j|1/r|k;l> ; |i> = MoCu(:,i).
!
!--- In this program , Two electron integral are express by using ------
!--- Chemist's notation------------------------------------------------- 
!-----------------------------------------------------------------------
    implicit none
    integer i,rec_char,narg,Node
    integer Mcyc,Conv,nstate,SCov,N_plot
    integer N,Nmem,Neri,Nocc,NoccA,NoccB,M,Me,Charge,Spin_Mult
    integer,allocatable:: bas_num(:),MAtom(:),NAtom(:),l_Basis(:,:),EI(:,:)
    real(kind=8) Ee,E_Nuc,mp_2,mp_3,Ecid,Ecisd
    real(kind=8) time_begin,time_end,time,tim(4)
    real(kind=8),allocatable:: En(:),MoCu(:,:),Density(:,:),Dipole(:,:,:)
    real(kind=8),allocatable:: EnA(:),MoCuA(:,:),DensityA(:,:)
    real(kind=8),allocatable:: EnB(:),MoCuB(:,:),DensityB(:,:)
    real(kind=8),allocatable:: Two_Electron(:,:),H_Core(:,:)
    real(kind=8),allocatable:: f_Mo(:,:),g_Mo(:,:)
    real(kind=8),allocatable:: f_MoA(:,:),f_MoB(:,:)
    real(kind=8),allocatable:: g_MoA(:,:),g_MoB(:,:),g_MoAB(:,:)
    real(kind=8),allocatable:: c_Basis(:,:,:),a_Basis(:,:,:),Rn(:,:),Rt(:,:)
    character(len=10) t1
    character(len=8) d1
    character(len=200) fname,sfname
    character(len=20) method,bas_set,option,outs,conver,maxcyc,A_or_B ! See ReadMe.txt
    character(len=20) diis,guess
    character(len=2),allocatable:: Atom(:),SAtom(:)
    character(len=4),allocatable:: IAtom(:),Ind_Basis(:)
    character(len=5) GridN
    logical exist1
    logical dft1,dft2,dft3,dft4,dft5,dft6,dft7,dft8,log1,log2
    
    call Head_Prient(t1,d1,1)
    !------------------ Read message from input file -------------------
    narg = IARGC()
    if( narg > 0 ) then
       call getarg(1,fname)
    else
       write(*,*) " input file name,or drag the file into terminal directly. "
       read(*,*) fname
    end if
    
    call Filename(fname,sfname,rec_char)
    open(unit = 10,file = fname)
    open(unit=100,file = sfname(1:rec_char)//".aop")
    call Head_Prient(t1,d1,2)
    read(10,*) M,Charge,Spin_Mult
    read(10,*) bas_set,method,option,outs
    read(10,*) maxcyc,Mcyc,conver,Conv,diis,guess
    read(10,*)
    write(*,*)
    write(*,*) " basis sets:  ",bas_set,";  method:  ",method
    write(*,*) " diis=",diis(6:8),"    ;    ","guess=",guess(7:12)
    !-------------------------------------------------------------------
    write(100,*) "------------------------- INPUT FILE ---------------------------"
    write(100,*)
    write(100,"(A2,G0,I4,I4)") "  ",M,Charge,Spin_Mult
    write(100,"(A2,A9,A12,A9,A6)") "  ",bas_set,method,option,outs
    write(100,"(A9,I4,A8,I3,A2,A12,A12)") "maxcyc=",Mcyc," conver=",Conv,"  ",diis,guess
    write(100,*)
    
    !-------------------------------------------------------------------
    allocate(Atom(M),MAtom(M),Rn(3,M))
             
    do i = 1,M
         read(10,*) Atom(i),Rn(1:3,i)
         write(100,"(A4,F20.8,2F14.8)") Atom(i),Rn(1:3,i)
    end do                                              !Angstroms
    write(100,*)
    write(100,*) "----------------------------------------------------------------" 
    write(100,*)
                                                    
    Rn  = Rn/0.5291772083d0                             !Bohr unit 
    !-------------------------------------------------------------------
    if( diis /= "diis=on" .and. diis /= "diis=off" ) then
        write(*,*) " Check diis option!"
        write(100,*) " Check diis option!"
        write(*,*)
        write(100,*)
        stop
    end if
    
    if( guess /= "guess=core" .and. guess /= "guess=atden") then
        write(*,*) " Check guess option!"
        write(100,*) " Check guess option!"
        write(*,*)
        write(100,*)
        stop
    end if
    !-------------------------------------------------------------------
    
    call Cal_NAtom(M,MAtom,Atom,Me)
    call Count_N(N,bas_set,MAtom,M)  
    !----------------------------------------------------------
    Nmem = 150 ! If N < Nmem ,store all Two_Electron in memory.
    if( N > Nmem ) then
        Neri = 1
    else
        Neri = N
    end if
    !----------------------------------------------------------
    write(*,*)
    write(*,"(2(A9,I5))") " Natom =",M," Nn_e =",Me
    write(100,"(2(A9,I5))") " Natom =",M," Nn_e =",Me

    allocate(Rt(3,N),EI(N,N),c_Basis(6,N,36),a_Basis(6,N,36),bas_num(N),NAtom(N),l_Basis(3,N))
    allocate(SAtom(N),IAtom(N),Ind_Basis(N),Dipole(N,N,3))
    
    !------------------- Charge and Spin Multiplicity ------------------
    if( abs(Charge) > M ) then
       write(*,*) " Segmentation Fault!"
       write(*,*) " Check Charge!"
       stop
    end if

    if ( mod(Spin_Mult,2) == 0 .and. mod(Me-Charge,2) == 0 ) then
       write(*,*) " Segmentation Fault!"
       write(*,*) " Check Charge and Spin Multiplicity!"
       stop
    end if
    
    if ( mod(Spin_Mult,2) == 1 .and. mod(Me-Charge,2) == 1 ) then
       write(*,*) " Segmentation Fault!"
       write(*,*) " Check Charge and Spin Multiplicity!"
       stop
    end if
        
    if ( Spin_Mult <= 0 ) then
       write(*,*) " Segmentation Fault!"
       write(*,*) " Check Charge and Spin Multiplicity!"
       stop
    end if

    if ( Spin_Mult == 1 ) then
        Nocc = (Me-Charge)/2
        write(*,"(3(A9,I5))") " Nocc  =",Nocc," Nbas =",N
        write(100,"(3(A9,I5))") " Nocc  =",Nocc," Nbas =",N
    end if
        
    if ( mod(Spin_Mult,2) == 0 ) then
        NoccA = (Me-Charge+Spin_Mult-1)/2
        NoccB = Me-Charge-NoccA
        write(*,"(3(A9,I5))") " NoccA =",NoccA," NoccB =",NoccB," Nbas =",N
        write(100,"(3(A9,I5))") " NoccA =",NoccA," NoccB =",NoccB," Nbas =",N
    end if
    
    if ( mod(Spin_Mult,2) == 1 .and. Spin_Mult > 1 ) then
        NoccA = (Me-Charge+Spin_Mult-1)/2
        NoccB = Me-Charge-NoccA
        write(*,"(3(A9,I5))") " NoccA =",NoccA," NoccB =",NoccB," Nbas =",N
        write(100,"(3(A9,I5))") " NoccA =",NoccA," NoccB =",NoccB," Nbas =",N
    end if
    write(*,*)
    !=========================== DFT Grid ==============================
    GridN = "grid2"
    
    if( GridN == "grid1" ) then    
        Node = 8760 !60*146
    end if
    
    if( GridN == "grid2" ) then 
        Node = 22650 !75*302
    end if
    !===================================================================
    write(100,*)
    write(100,*)
    call Cpu_time(time_begin)
    
    if ( Spin_Mult == 1 .and. mod(Me-Charge,2) == 0 ) then
    
       !================== Close Shell,Multiplicity = 1 ===============!
       
       !----------------------------- HF -------------------------------
       allocate(En(N),MoCu(N,N),Density(N,N),Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2), &
                H_Core(N,N))
                
       log1 = method == "hf" .and. option == "none"
       log2 = method == "rhf" .and. option == "none"
       
       if (log1 .or. log2) then
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
       end if
       !----------------------------- DFT ------------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       if ( log1 .and. option == "none" ) then
          call DFT_Close(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                         Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                         c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
       end if
       !-------------------------- HF mobt -----------------------------
       log1 = method == "hf" .and. option == "mobt"
       log2 = method == "rhf" .and. option == "mobt"
       
       if (log1 .or. log2) then
       
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) En
              read(555,*) 
              read(555,*) MoCu
              
          else
              call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                       Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                       c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)         
              write(555,*) En
              write(555,*) 
              write(555,*) MoCu
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
                  
          write(*,*) " Press d to draw Orbital:"
          read(*,*) A_or_B
          if ( A_or_B == "d" ) call Molecular_orbital(N,MoCu,En,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !-------------------------- DFT mobt ----------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       
       if ( log1 .and. option == "mobt" ) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) En
              read(555,*) 
              read(555,*) MoCu
              
          else
              call DFT_Close(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                             Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                             c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)         
              write(555,*) En
              write(555,*) 
              write(555,*) MoCu
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if

          write(*,*) " Press d to draw Orbital:"
          read(*,*) A_or_B
          if ( A_or_B == "d" ) call Molecular_orbital(N,MoCu,En,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !-------------------------- HF den ------------------------------
       log1 = method == "hf" .and. option == "den"
       log2 = method == "rhf" .and. option == "den"
       
       if (log1 .or. log2) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) En
              read(555,*) 
              read(555,*) MoCu
              Density = 2.0d0*matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
          else
              call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                       Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                       c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)         
              write(555,*) En
              write(555,*) 
              write(555,*) MoCu
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
          
          call Electron_Density(N,Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !-------------------------- DFT den -----------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       
       if ( log1 .and. option == "den" ) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) En
              read(555,*) 
              read(555,*) MoCu
              Density = 2.0d0*matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
          else
              call DFT_Close(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                             Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                             c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)         
              write(555,*) En
              write(555,*) 
              write(555,*) MoCu
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
          
          call Electron_Density(N,Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !---------------------------- MP2 -------------------------------
       log1 = method == "mp2" .and. option == "none"
       log2 = method == "rmp2" .and. option == "none"
       
       if ( log1 .or. log2 ) then
          allocate(f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2))
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,"N")
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCu,Density,Two_Electron,H_Core,f_Mo)
          
          call MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
       end if
       !---------------------------- MP3 -------------------------------
       log1 = method == "mp3" .and. option == "none"
       log2 = method == "rmp3" .and. option == "none"
       
       if ( log1 .or. log2 ) then
          allocate(f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2))
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,"N")
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCu,Density,Two_Electron,H_Core,f_Mo)
          
          call MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
          call MP3_Close(N,Nocc,En,g_Mo,EI,mp_3)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(*,"(A9,F14.7)") "MP3  = ",Ee + E_Nuc  + mp_2 + mp_3
          write(100,"(A9,F14.7)") "MP3  = ",Ee + E_Nuc + mp_2 + mp_3
       end if
       !---------------------------- TDHF ------------------------------
       if (method == "tdhf" ) then
                 
          do i = 1,10
              if( option(i:i) == "#" ) exit
          end do
          if( i > 9 ) then
              write(*,*) " Abnormal termination, please check the option!"
              write(100,*) " Abnormal termination, please check the option!"
              write(*,*)
              write(100,*)
              stop
          end if
          
          allocate(f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2))
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,"N")
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(f_Mo,Density,Two_Electron,H_Core)
          
          call MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
          
          read( option(1:i-1),"(I10)" ) nstate
          read( option(i+1:10),"(I10)" ) SCov
          
          call TDHF_Close(nstate,SCov,En,Ee,g_Mo,EI,MoCu,Dipole,N,Nocc,E_Nuc)
          
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
       end if
       !---------------------------- CIS -------------------------------
       if (method == "cis" ) then
                 
          do i = 1,10
              if( option(i:i) == "#" ) exit
          end do
          if( i > 9 ) then
              write(*,*) " Abnormal termination, please check the option!"
              write(100,*) " Abnormal termination, please check the option!"
              write(*,*)
              write(100,*)
              stop
          end if
          
          allocate(f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2))
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,"N")
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(f_Mo,Density,Two_Electron,H_Core)
          
          call MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
          
          read( option(1:i-1),"(I10)" ) nstate
          read( option(i+1:10),"(I10)" ) SCov
          
          call CIS_Close(nstate,SCov,En,Ee,g_Mo,EI,MoCu,Dipole,N,Nocc,E_Nuc)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
       end if
       !---------------------------- CID -------------------------------
       if (method == "cid".and. option == "none") then
          allocate(f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2))
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,"W")
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCu,Density,Two_Electron,H_Core)
          
          call MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
          call CID_Close(f_Mo,g_Mo,EI,N,Nocc,E_Nuc,Ecid)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(*,"(A9,F14.7)") "CID  = ",Ecid + E_Nuc
          write(100,"(A9,F14.7)") "CID  = ",Ecid + E_Nuc
       end if
       !---------------------------- CISD ------------------------------
       if (method == "cisd".and. option == "none") then
          allocate(f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2))
          call RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                   c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,"W")
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCu,Density,Two_Electron,H_Core)
          
          call MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
          call CISD_Close(f_Mo,g_Mo,EI,N,Nocc,E_Nuc,Ecisd)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(*,"(A9,F14.7)") "CISD = ",Ecisd + E_Nuc
          write(100,"(A9,F14.7)") "CISD = ",Ecisd + E_Nuc
       end if
    !===================================================================
    !===================================================================
    else
       !=================== Open Shell,Multiplicity > 1 ===============!
       
       !----------------------------- HF -------------------------------
       allocate(EnA(N),MoCuA(N,N),DensityA(N,N),EnB(N),MoCuB(N,N),DensityB(N,N), &
                Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2),H_Core(N,N))
                
       log1 = method == "hf" .and. option == "none"
       log2 = method == "uhf" .and. option == "none"
       
       if (log1 .or. log2) then
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
       end if 
       !----------------------------- ROHF -----------------------------
       if (method == "rohf".and. option == "none") then
          deallocate(MoCuA,MoCuB,EnA,EnB)
          allocate(MoCu(N,N),En(N))
          call ROHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                    NAtom,Two_Electron,EI,H_Core,MoCu,Dipole,&
                    En,Ee,DensityA,DensityB,c_Basis,a_Basis,&
                    l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
       end if 
       !----------------------------- DFT ------------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       if ( log1 .and. option == "none" ) then
          call DFT_Open(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                        NAtom,Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                        EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                        l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
       end if
       !---------------------------- RODFT -----------------------------
       dft1 = method == "roxalpha" 
       dft2 = method == "ropbepbe" 
       dft3 = method == "rob3lyp" 
       dft4 = method == "roblyp" 
       dft5 = method == "ropbe1pbe" 
       dft6 = method == "rolsda" 
       dft7 = method == "rob3pw91" 
       dft8 = method == "ropw91pw91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       if ( log1 .and. option == "none" ) then
          deallocate(MoCuA,MoCuB,EnA,EnB)
          allocate(MoCu(N,N),En(N))
          call DFT_RO(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                    NAtom,Two_Electron,EI,H_Core,MoCu,Dipole,&
                    En,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                    l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
       end if
       !-------------------------- HF mobt -----------------------------
       log1 = method == "hf" .and. option == "mobt"
       log2 = method == "uhf" .and. option == "mobt"
       
       if (log1 .or. log2) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) EnA
              read(555,*)
              read(555,*) Enb
              read(555,*)  
              read(555,*) MoCuA
              read(555,*)  
              read(555,*) MoCuB
                            
          else
              call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                       Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                       EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                       l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)                    
              write(555,*) EnA
              write(555,*)
              write(555,*) Enb
              write(555,*)  
              write(555,*) MoCuA
              write(555,*)  
              write(555,*) MoCuB
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if

          write(*,*)
          write(*,*) " Press a to draw Alpha Orbital,or b to draw Beta Orbital:"
          read(*,*) A_or_B
          if ( A_or_B == "a" ) call Molecular_Orbital(N,MoCuA,EnA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
          if ( A_or_B == "b" ) call Molecular_Orbital(N,MoCuB,EnB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)        
       end if
       !-------------------------- DFT mobt ----------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       
       if ( log1 .and. option == "mobt" ) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) EnA
              read(555,*)
              read(555,*) Enb
              read(555,*)  
              read(555,*) MoCuA
              read(555,*)  
              read(555,*) MoCuB
                            
          else
              call DFT_Open(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                            NAtom,Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                            EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                            l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)                    
              write(555,*) EnA
              write(555,*)
              write(555,*) Enb
              write(555,*)  
              write(555,*) MoCuA
              write(555,*)  
              write(555,*) MoCuB
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
          
          write(*,*)
          write(*,*) " Press a to draw Alpha Orbital,or b to draw Beta Orbital:"
          read(*,*) A_or_B
          if ( A_or_B == "a" ) call Molecular_Orbital(N,MoCuA,EnA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
          if ( A_or_B == "b" ) call Molecular_Orbital(N,MoCuB,EnB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)   
       end if
       !-------------------------- HF den ------------------------------
       log1 = method == "hf" .and. option == "den"
       log2 = method == "uhf" .and. option == "den"
       
       if (log1 .or. log2) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) EnA
              read(555,*)
              read(555,*) Enb
              read(555,*)  
              read(555,*) MoCuA
              read(555,*)  
              read(555,*) MoCuB
              DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
              DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB)))                           
          else
              call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                       Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                       EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                       l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)                    
              write(555,*) EnA
              write(555,*)
              write(555,*) Enb
              write(555,*)  
              write(555,*) MoCuA
              write(555,*)  
              write(555,*) MoCuB
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
          
          call Electron_Density(N,DensityA+DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !-------------------------- DFT den -----------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       
       if ( log1 .and. option == "den" ) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) EnA
              read(555,*)
              read(555,*) Enb
              read(555,*)  
              read(555,*) MoCuA
              read(555,*)  
              read(555,*) MoCuB
              DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
              DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB)))                           
          else
              call DFT_Open(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                            NAtom,Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                            EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                            l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)                    
              write(555,*) EnA
              write(555,*)
              write(555,*) Enb
              write(555,*)  
              write(555,*) MoCuA
              write(555,*)  
              write(555,*) MoCuB
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
          
          call Electron_Density(N,DensityA+DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num) 
       end if
       !-------------------------- HF sden -----------------------------
       log1 = method == "hf" .and. option == "sden"
       log2 = method == "uhf" .and. option == "sden"
       
       if (log1 .or. log2) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) EnA
              read(555,*)
              read(555,*) Enb
              read(555,*)  
              read(555,*) MoCuA
              read(555,*)  
              read(555,*) MoCuB
              DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
              DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB)))                           
          else
              call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                       Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                       EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                       l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)                    
              write(555,*) EnA
              write(555,*)
              write(555,*) Enb
              write(555,*)  
              write(555,*) MoCuA
              write(555,*)  
              write(555,*) MoCuB
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if
          call Electron_Density(N,DensityA-DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !-------------------------- DFT sden ----------------------------
       dft1 = method == "XAlpha" .or. method == "xalpha"
       dft2 = method == "PBEPBE" .or. method == "pbepbe"
       dft3 = method == "b3lyp" .or. method == "B3LYP"
       dft4 = method == "blyp" .or. method == "BLYP"
       dft5 = method == "PBE1PBE" .or. method == "pbe1pbe"
       dft6 = method == "LSDA" .or. method == "lsda"
       dft7 = method == "b3pw91" .or. method == "B3PW91"
       dft8 = method == "pw91pw91" .or. method == "PW91PW91"
       log1 = dft1 .or. dft2 .or. dft3 .or. dft4 .or. dft5 .or. dft6 .or. dft7 .or. dft8
       
       if ( log1 .and. option == "sden" ) then
          inquire(exist=exist1,file=sfname(1:rec_char)//".mo")
          open(unit=555,file=sfname(1:rec_char)//".mo")
          
          if( exist1 ) then
              call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                   bas_num,Atom,IAtom,SAtom,Ind_Basis)
              write(*,*) " Using molecular orbital file."    
              write(*,*)                   
              read(555,*) N_plot
              if( N_plot /= N ) stop
              read(555,*)                    
              read(555,*) EnA
              read(555,*)
              read(555,*) Enb
              read(555,*)  
              read(555,*) MoCuA
              read(555,*)  
              read(555,*) MoCuB
              DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
              DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB)))                           
          else
              call DFT_Open(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                            NAtom,Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                            EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                            l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
                       
              write(555,*) N
              write(555,*)                    
              write(555,*) EnA
              write(555,*)
              write(555,*) Enb
              write(555,*)  
              write(555,*) MoCuA
              write(555,*)  
              write(555,*) MoCuB
              
              write(*,*)
              write(100,*)
              write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
              write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          end if

          call Electron_Density(N,DensityA-DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
       end if
       !---------------------------- MP2 -------------------------------     
       log1 = method == "mp2" .and. option == "none"
       log2 = method == "ump2" .and. option == "none"
       
       if ( log1 .or. log2 ) then
          allocate(g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N))
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuA,f_Mo,g_MoA,"N")
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuB,f_Mo,g_MoB,"N")
          call Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_MoAB)                
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCuA,DensityA,MoCuB,DensityB,H_Core,Two_Electron,f_MoA,f_MoB)
          
          call MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
       end if
       !---------------------------- MP3 ------------------------------- 
       log1 = method == "mp3" .and. option == "none"
       log2 = method == "ump3" .and. option == "none"
       
       if ( log1 .or. log2 ) then
          allocate(g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N))
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuA,f_Mo,g_MoA,"N")
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuB,f_Mo,g_MoB,"N")
          call Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_MoAB)
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCuA,DensityA,MoCuB,DensityB,H_Core,Two_Electron,f_MoA,f_MoB)
          
          call MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
          call MP3_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_3)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(*,"(A9,F14.7)") "MP3  = ",Ee + E_Nuc  + mp_2 + mp_3
          write(100,"(A9,F14.7)") "MP3  = ",Ee + E_Nuc + mp_2 + mp_3
       end if
       !---------------------------- TDHF ------------------------------ 
       if (method == "tdhf") then
       
          do i = 1,10
              if( option(i:i) == "#" ) exit
          end do
          if( i > 9 ) then
              write(*,*) " Abnormal termination, please check the option!"
              write(100,*) " Abnormal termination, please check the option!"
              write(*,*)
              write(100,*)
              stop
          end if
          
          allocate(g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N))
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuA,f_MoA,g_MoA,"N")
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuB,f_MoB,g_MoB,"N")
          call Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_MoAB)
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(f_MoA,f_MoB,DensityA,DensityB,H_Core,Two_Electron)
          
          call MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
          
          read( option(1:i-1),"(I10)" ) nstate
          read( option(i+1:10),"(I10)" ) SCov
          
          call TDHF_Open(nstate,SCov,EnA,EnB,Ee,g_MoA,g_MoB,g_MoAB,EI,MoCuA, &
                         MoCuB,Dipole,N,NoccA,NoccB,E_Nuc)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
       end if
       !---------------------------- CIS ------------------------------- 
       if (method == "cis") then
       
          do i = 1,10
              if( option(i:i) == "#" ) exit
          end do
          if( i > 9 ) then
              write(*,*) " Abnormal termination, please check the option!"
              write(100,*) " Abnormal termination, please check the option!"
              write(*,*)
              write(100,*)
              stop
          end if
          
          allocate(g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N))
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuA,f_MoA,g_MoA,"N")
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuB,f_MoB,g_MoB,"N")
          call Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_MoAB)
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(f_MoA,f_MoB,DensityA,DensityB,H_Core,Two_Electron)
          
          call MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
          
          read( option(1:i-1),"(I10)" ) nstate
          read( option(i+1:10),"(I10)" ) SCov
          
          call CIS_Open(nstate,SCov,EnA,EnB,Ee,g_MoA,g_MoB,g_MoAB,EI,MoCuA, &
                        MoCuB,Dipole,N,NoccA,NoccB,E_Nuc)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
       end if
       !---------------------------- CID -------------------------------  
       if (method == "cid" .and. option == "none") then
          allocate(g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N))
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuA,f_MoA,g_MoA,"Y")
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuB,f_MoB,g_MoB,"W")
          call Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_MoAB)
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCuA,DensityA,MoCuB,DensityB,H_Core,Two_Electron)
          
          call MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
          call CID_Open(f_MoA,f_MoB,g_MoA,g_MoB,g_MoAB,EI,N,NoccA,NoccB,E_Nuc,Ecid)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(*,"(A9,F14.7)") "CID  = ",Ecid + E_Nuc
          write(100,"(A9,F14.7)") "CID  = ",Ecid + E_Nuc
       end if
       !---------------------------- CISD ------------------------------ 
       if (method == "cisd" .and. option == "none") then
          allocate(g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N))
          call UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom,NAtom, &
                   Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                   EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                   l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
                   
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuA,f_MoA,g_MoA,"Y")
          call Mol_Obt(N,Two_Electron,EI,H_Core,MoCuB,f_MoB,g_MoB,"W")
          call Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_MoAB)
          write(*,*) " <ia|jb> transformation complete."
          
          deallocate(MoCuA,DensityA,MoCuB,DensityB,H_Core,Two_Electron)
          
          call MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
          call CISD_Open(f_MoA,f_MoB,g_MoA,g_MoB,g_MoAB,EI,N,NoccA,NoccB,E_Nuc,Ecisd)
          write(*,*)
          write(100,*)
          write(*,"(A9,F14.7)") "HF   = ",Ee + E_Nuc 
          write(100,"(A9,F14.7)") "HF   = ",Ee + E_Nuc
          write(*,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(100,"(A9,F14.7)") "MP2  = ",Ee + E_Nuc  + mp_2
          write(*,"(A9,F14.7)") "CISD = ",Ecisd + E_Nuc
          write(100,"(A9,F14.7)") "CISD = ",Ecisd + E_Nuc
       end if
       !----------------------------------------------------------------
       
    end if
    
    call Cpu_time(time_end)
    time = time_end - time_begin
    
    write(*,*) 
    write(100,*)
    call Cotime(time,tim)
    call Tail_Prient(t1,d1) 
    call Time_Prient(tim)
    
    if ( time > 0.0001d0 ) then
       write(*,*) "Normal termination, the program has finished now!"
       write(100,*) "Normal termination, the program has finished now!"
       write(*,*)
       write(100,*)
    else
       write(*,*) "Abnormal termination, please check the method and option!"
       write(100,*) "Abnormal termination, please check the method and option!"
       write(*,*)
       write(100,*)
    end if   
        
    write(100,*) "===================================================================="    
    close(100)
    
end subroutine Amesp
!=======================================================================
subroutine Cal_NAtom(M,MAtom,Atom,Me)
    ! MAtom,Atom and Me.
    implicit none
    integer M,MAtom(M),i,j,Me
    character(len=2) Atom(M),Atom_Set(103)
    
    Atom_Set = &
  (/'H ',                                                                                'He', &  !   2
    'Li','Be',                                                  'B ','C ','N ','O ','F ','Ne', &  !  10
    'Na','Mg',                                                  'Al','Si','P ','S ','Cl','Ar', &  !  18
    'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr', &  !  36
    'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe', &  !  54
    'Cs','Ba',                                                                                 &  !  56
    'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',                     &  !  70
              'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn', &  !  86
    'Fr','Ra',                                                                                 &  !  88
              'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',           &  ! 102
    'Lr'                                                                                       &  ! 103
  /)
    
    do i = 1,M
    
       do j = 1,36
          if( Atom(i) == Atom_Set(j) ) exit
       end do
       MAtom(i) = j
       
       if( j > 36 ) then
           write(*,*) " Warning : Only support H--Kr now!"
           stop
       end if
       
    end do
    
    Me = sum(MAtom)
        
end subroutine Cal_NAtom
!=======================================================================
subroutine Count_N(N,bas_set,MAtom,M)
    ! Calculate the number of basis.
    implicit none
    integer i,M,N,MAtom(M),N_bas(36)
    character(len=10) bas_set
    
    if( bas_set == "sto-3g" ) then
    
        N_bas = (/1,1,5,5,5,5,5,5,5,5,9,9,9,9,9,9,9,9,13,13,19,19,19,19,19, &
                 19,19,19,19,19,19,19,19,19,19,19/)
        N = 0
        do i = 1,M
           N = N + N_bas(MAtom(i))
        end do

    else if( bas_set == "3-21g"  ) then
    
        N_bas = (/2,2,9,9,9,9,9,9,9,9,13,13,13,13,13,13,13,13,17,17,29,29, &
                 29,29,29,29,29,29,29,29,23,23,23,23,23,23/)
        N = 0
        do i = 1,M
           N = N + N_bas(MAtom(i))
        end do
        
    else if( bas_set == "6-31g"  ) then
    
        N_bas = (/2,2,9,9,9,9,9,9,9,9,13,13,13,13,13,13,13,13,17,17,29,29, &
                 29,29,29,29,29,29,29,29,24,24,24,24,24,24/)
        N = 0
        do i = 1,M
           N = N + N_bas(MAtom(i))
        end do
        
    else if( bas_set == "6-31g*"  ) then
    
        N_bas(1:20) = (/2,2,15,15,15,15,15,15,15,15,19,19,19, &
                       19,19,19,19,19,23,23/)
        N = 0
        do i = 1,M
           N = N + N_bas(MAtom(i))
        end do
        
    else if( bas_set == "6-31g**"  ) then
    
        N_bas(1:20) = (/5,5,15,15,15,15,15,15,15,15,19,19,19, &
                       19,19,19,19,19,23,23/)
        N = 0
        do i = 1,M
           N = N + N_bas(MAtom(i))
        end do
        
    else if( bas_set == "6-311g"  ) then
    
        N_bas(1:20) = (/3,3,13,13,13,13,13,13,13,13,21,21,21,21,21,21,21,21,35,35/)
        
        N = 0
        do i = 1,M
           N = N + N_bas(MAtom(i))
        end do
        
    else
       write(100,*)
       write(*,*) " Abnormal termination, please check the Basis Sets!"
       write(100,*) " Abnormal termination, please check the Basis Sets!"
       write(*,*)
       write(100,*)
       stop
    end if
    
end subroutine Count_N
!=======================================================================

