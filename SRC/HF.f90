subroutine RHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
              Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
              c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
!------------------Close Shell，Multiplicity = 1------------------------

    implicit none
    integer Mcyc,Conv
    integer N,Nmem,Neri,Nocc,M,bas_num(N),l_Basis(3,N)
    integer i,j,ii,idi,MAtom(M),NAtom(N),EI(N,N)
    real(kind=8) Ee,E_Nuc,E_T,E_V,Delta_E,Delta_Den
    real(kind=8) Rn(3,M),Rt(3,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) MoCu(N,N),En(N),Density(N,N),Mulliken(N,N)
    real(kind=8) Dipole(N,N,3),Overlap(N,N),kinetic(N,N),Density_T(N,N)
    real(kind=8) N_E_Attract(N,N),H_Core(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) Fock(N,N,0:9),G(N,N),Sqrt_S(N,N)
    real(kind=8) Er(N,N,9),emax,CCoff(10),ethre
    real(kind=8) Temp1,Temp2,ETemp,Overlap_Temp(N,N)
    real(kind=8) dipex,dipey,dipez,dipnx,dipny,dipnz,dipx,dipy,dipz,diptot
    character(len=1) Sn(N)
    character(len=20) outs,bas_set,diis,guess
    character(len=2) Atom(M),SAtom(N)
    character(len=4) Ind_Basis(N),IAtom(N)
    character(len=10) t2
    character(len=8) d2
    logical jud_diis
    
    Sn(1:Nocc) = "O" ; Sn(Nocc+1:N) = "V"   
    jud_diis = .true.

    call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                    bas_num,Atom,IAtom,SAtom,Ind_Basis)
    
    !-----------------------One_Electron_Int----------------------------
    call Overlap_Kinetic_Int(Overlap,Kinetic,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)
    call N_E_Attractive_Int(N_E_Attract,M,N,a_Basis,c_Basis,l_Basis,bas_num,Rn,Rt,MAtom,NAtom)
    call Dipole_Int(Dipole,Overlap,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)   
     
    write(*,*) " One-electron integral complete."
     
    if( outs == "out=3" .or. outs == "out=4" ) then  
       write(100,*) "===== Overlap Matrix: ====="
       call Mat_Out1(Overlap,N)
       write(100,*)
       write(100,*)
    
       write(100,*) "===== Kinetic Matrix: ====="
       call Mat_Out1(Kinetic,N)
       write(100,*)
       write(100,*)

       write(100,*) "===== Nuclear Attraction: ====="
       call Mat_Out1(N_E_Attract,N)
       write(100,*)
       write(100,*)
       
       write(100,*) "===== Multipole matrices(X): ====="
       call Mat_Out1(Dipole(:,:,1),N)
       write(100,*)
       write(100,*)
       
       write(100,*) "===== Multipole matrices(Y): ====="
       call Mat_Out1(Dipole(:,:,2),N)
       write(100,*)
       write(100,*)
             
       write(100,*) "===== Multipole matrices(Z): ====="
       call Mat_Out1(Dipole(:,:,3),N)
       write(100,*)
       write(100,*)
    end if
    !-----------------------Two_Electron_Int----------------------------
    if ( N > Nmem ) then
    
        call Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)
         
        write(*,*) " Two-electron integral complete."
        write(*,*)  
        
        if( outs == "out=4" ) then
            write(100,*) " Too many Two-electron integral !"
            write(100,*)
        end if
          
    else
    
        call Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)
        write(*,*) " Two-electron integral complete."
        write(*,*)
        
        if( outs == "out=4" ) then  
            call Two_Ele_Int_Out(N,Neri,Two_Electron,EI)
        end if
        
    end if
    !------------------------Energy between Nucleus---------------------    
    
    call Nucleus_Energy(M,Rn,MAtom,E_Nuc)
    
    !-------------------------SCF Procedure-----------------------------
    write(*,*) "=================== SCF Procedure =====================" 
    write(100,*) "=================== SCF Procedure ====================="
    write(*,*)
    write(100,*)
    ETemp = 10000000.0d0

    H_Core =  Kinetic + N_E_Attract
    Overlap_Temp = Overlap
    call Mat_Sqrt(Overlap_Temp,Sqrt_S,N)
    
    !-------------------------------------------------------------------
    if (guess == "guess=core" ) then
         call Main_Calculate_S(Sqrt_S,H_Core,MoCu,En,N) !Initial Molecular Orbit.
    end if
    
    if (guess == "guess=atden" ) then
         call AtDen(M,N,Nmem,Neri,MAtom,NAtom,Overlap,Kinetic,H_core,Rn,Rt,a_Basis,c_Basis, &
                    l_Basis,bas_num,Two_Electron,EI,bas_set,t2,d2,MoCu) !Initial Molecular Orbit.
    end if
    write(*,*) "                  INITIAL GUESS DONE"
    write(*,*) 
    !-------------------------------------------------------------------
    Density = matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
    
    do ii = 1,Mcyc                !The maximum number of cycles
       write(*,"(A7,I4)") "Cycle",ii
       write(100,"(A7,I4)") "Cycle",ii
       
       Density_T = Density
       Density = matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
       
       Temp1 = 0.0d0
       do i = 1,N 
          Temp2 = 0.0d0
          do j = 1,N
             Temp2 = Temp2 + (Density_T(i,j)-Density(i,j))**2
          end do
          Temp1 = Temp1 + Temp2
       end do
       
       Delta_Den = dsqrt(Temp1)/real(N,8)
          
       if(jud_diis) then
           if( ii < 5 ) Density = 0.2d0*Density_T + 0.8d0*Density
           if( ii >=5 ) Density = 0.5d0*Density_T + 0.5d0*Density
       end if
       
       !------------------------ G matrix ------------------------------
       call RHF_Const_G(N,Nmem,Neri,Two_Electron,EI,Density,G,t2,d2)
       !----------------------------------------------------------------
    
       idi = mod(ii,9)
       if( mod(ii,9) == 0 ) idi = 9
       
       Fock(:,:,idi) = H_core + G         !Fock Matrix
       
       if ( diis == "diis=on" ) then 
       !---------------------------  DIIS  -----------------------------
       call Error_Mat(N,Fock(:,:,idi),Overlap,Density,Er(:,:,idi),emax)
       
       if( ii < 5 .and. jud_diis ) ethre = 1.0d-1
       if( ii >=5 .and. jud_diis ) ethre = 5.0d-1
       
       if( emax < ethre) then
       
           jud_diis = .false.
           if( ii <= 9 ) then
                call C_Coff(idi,N,Er(:,:,1:idi),CCoff(1:idi+1))
                Fock(:,:,0) = 0.0d0
                do i = 1,idi
                   Fock(:,:,0) = Fock(:,:,0) + Fock(:,:,i)*CCoff(i)
                end do 
           else
                call C_Coff(9,N,Er(:,:,1:9),CCoff(1:10))
                Fock(:,:,0) = 0.0d0
                do i = 1,9
                   Fock(:,:,0) = Fock(:,:,0) + Fock(:,:,i)*CCoff(i)
                end do 
           end if
           
       else
           Fock(:,:,0) =  Fock(:,:,idi)    
       end if
       !----------------------------------------------------------------
       else
          emax = 0.0d0
          Fock(:,:,0) =  Fock(:,:,idi)
       end if
       
       call Main_Calculate_S(Sqrt_S,Fock(:,:,0),MoCu,En,N) !Solve FC=SCE

       Temp1 = 0.0d0
       do i = 1,N
           Temp2 = 0.0d0
           do j = 1,N
               Temp2 = Temp2 + Density(j,i)*(H_core(j,i) + Fock(j,i,0))
           end do
           Temp1 = Temp1 + Temp2 
       end do

       Ee = Temp1                 !Calculate Energy E
       Delta_E = dabs( Ee - ETemp )
       
       write(*,"(A14,ES16.8)") "E(HF) = ",Ee + E_Nuc 
       write(*,"(A14,ES14.6)") "DIIS_ER=",emax
       write(*,"(A13,ES15.6,A10,ES16.6)") "DDen  =",Delta_Den,"DE =",Delta_E
       write(100,"(A14,ES16.8)") "E(HF) = ",Ee + E_Nuc 
       write(100,"(A14,ES14.6)") "DIIS_ER=",emax
       write(100,"(A13,ES15.6,A10,ES16.6)") "DDen  =",Delta_Den,"DE =",Delta_E
       write(*,*)
       write(100,*)
       if( Delta_E < 10.0d0**(-Conv) .and. Delta_Den < 1.0d-4 ) exit 
       !Judge whether the SCF procedure is converged.
       ETemp = Ee

    end do
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do i = 1,N
        Temp2 = 0.0d0
        do j = 1,N
            Temp2 = Temp2 + 2.0d0*Density(j,i)*Kinetic(j,i)
        end do
        Temp1 = Temp1 + Temp2 
    end do  
    E_T = Temp1 ; E_V = Ee + E_Nuc - E_T
    !-------------------------------------------------------------------  
    write(*,*)"======================================================="
    write(100,*)"======================================================="

    if( ii <= Mcyc-1 ) then
        write(*,"(A18,I4,A9)") " SCF Done , After",ii," cycles."
        write(100,"(A18,I4,A9)") " SCF Done , After",ii," cycles."
        write(100,"(A9,F7.4,A9,ES14.6)") "  -V/T = ",-E_V/E_T,"    KE =",E_T
        write(*,*)"======================================================="
       write(100,*)"======================================================="
    end if 
    
    if( ii > Mcyc-1 ) then
        write(100,*)
        write(100,*) " Fail to convergence!"
        write(*,*)
        write(*,*) " Fail to convergence!"
        stop
    end if

    !----------------------------Print----------------------------------
    
    if( outs == "out=1" .or. outs == "out=2" .or. outs == "out=3" .or. outs == "out=4") then
       write(*,*)
       write(100,*)
       write(100,*) "Molecular Orbital Coefficients:"
       write(100,*)    
       call Mat_Out2(MoCu,En,Sn,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
    end if
    
    if( outs == "out=2" .or. outs == "out=3" .or. outs == "out=4") then 
    
      Density = matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
      Density = 2.0d0*Density   !Density Matrix 
      
      do i = 1,N
          do j = i,N
             Temp1 = Density(j,i)*Overlap(j,i)
             Mulliken(i,j) = Temp1
             Mulliken(j,i) = Temp1
          end do
       end do                    !Full Mulliken population analysis
       
       !------------------------- Dipole Moment ------------------------  
       dipex = sum(Density*Dipole(:,:,1))
       dipey = sum(Density*Dipole(:,:,2))
       dipez = sum(Density*Dipole(:,:,3))
       
       dipnx = 0.0d0 ; dipny = 0.0d0 ; dipnz = 0.0d0
       do i = 1,M
           dipnx = dipnx + MAtom(i)*Rn(1,i)
           dipny = dipny + MAtom(i)*Rn(2,i)
           dipnz = dipnz + MAtom(i)*Rn(3,i)
       end do
       
       dipx = (- dipex + dipnx)*2.54177d0
       dipy = (- dipey + dipny)*2.54177d0
       dipz = (- dipez + dipnz)*2.54177d0
       diptot = dsqrt( dipx**2 + dipy**2 + dipz**2 )
       !---------------------------------------------------------------- 
          
       write(100,*) "Density Matrix:"
       call Mat_Out1(Density,N)
       write(100,*)
       
       write(100,*) "Full Mulliken population analysis:"
       call Mat_Out1(Mulliken,N)
       write(100,*)
       
       write(100,*) "Gross orbital populations:"
       write(100,*)
       call GOP_Out1(Mulliken,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
       
       write(100,*) "Dipole moment (field-independent basis, Debye):"
       write(100,*)
       write(100,"(3(A6,F10.4),A8,F10.4)") "X=",dipx,"Y=",dipy,"Z=",dipz,"Tot=",diptot
       write(100,*)
       
    end if
    
end subroutine RHF
!=======================================================================
!=======================================================================
subroutine UHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
               NAtom,Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
               EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
               l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
!------------------Open Shell，Multiplicity >= 2------------------------

    implicit none
    integer Mcyc,Conv,idi
    integer N,Nmem,Neri,NoccA,NoccB,M,bas_num(N),l_Basis(3,N)
    integer i,j,ii,MAtom(M),NAtom(N),EI(N,N)
    real(kind=8) Ee,E_Nuc,E_T,E_V,Delta_E,S2
    real(kind=8) Rn(3,M),Rt(3,N),Coul(N,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) MoCuA(N,N),EnA(N),DensityA(N,N),MullikenA(N,N)
    real(kind=8) MoCuB(N,N),EnB(N),DensityB(N,N),MullikenB(N,N)
    real(kind=8) DensityA_T(N,N),DensityB_T(N,N)
    real(kind=8) Dipole(N,N,3),Overlap(N,N),kinetic(N,N)
    real(kind=8) N_E_Attract(N,N),H_Core(N,N),Sqrt_S(N,N)
    real(kind=8) Overlap_Temp(N,N),Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) dipex,dipey,dipez,dipnx,dipny,dipnz,dipx,dipy,dipz,diptot
    real(kind=8) FockA(N,N,0:9),ExcA(N,N),FockB(N,N,0:9),ExcB(N,N)
    real(kind=8) ErA(N,N,9),emaxA,CCoff(10)
    real(kind=8) ErB(N,N,9),emaxB,ethre
    real(kind=8) Temp1,Temp2,EeA,EeB,ETemp,Delta_DenA,Delta_DenB
    character(len=1) SnA(N),SnB(N)
    character(len=20) outs,bas_set
    character(len=20) diis,guess
    character(len=2) Atom(M),SAtom(N)
    character(len=4) Ind_Basis(N)
    character(len=4) IAtom(N)
    character(len=10) t2
    character(len=8) d2
    logical jud_diis
   
    SnA(1:NoccA) = "O" ; SnA(NoccA+1:N) = "V" 
    SnB(1:NoccB) = "O" ; SnB(NoccB+1:N) = "V" 
    jud_diis = .true.
    
    call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                    bas_num,Atom,IAtom,SAtom,Ind_Basis)
    
    !-----------------------One_Electron_Int----------------------------
    call Overlap_Kinetic_Int(Overlap,Kinetic,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)
    call N_E_Attractive_Int(N_E_Attract,M,N,a_Basis,c_Basis,l_Basis,bas_num,Rn,Rt,MAtom,NAtom)
    call Dipole_Int(Dipole,Overlap,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)   
     
    write(*,*) " One-electron integral complete."
     
    if( outs == "out=3" .or. outs == "out=4" ) then  
       write(100,*) "===== Overlap Matrix: ====="
       call Mat_Out1(Overlap,N)
       write(100,*)
       write(100,*)
    
       write(100,*) "===== Kinetic Matrix: ====="
       call Mat_Out1(Kinetic,N)
       write(100,*)
       write(100,*)

       write(100,*) "===== Nuclear Attraction: ====="
       call Mat_Out1(N_E_Attract,N)
       write(100,*)
       write(100,*)
       
       write(100,*) "===== Multipole matrices(X): ====="
       call Mat_Out1(Dipole(:,:,1),N)
       write(100,*)
       write(100,*)
       
       write(100,*) "===== Multipole matrices(Y): ====="
       call Mat_Out1(Dipole(:,:,2),N)
       write(100,*)
       write(100,*)
             
       write(100,*) "===== Multipole matrices(Z): ====="
       call Mat_Out1(Dipole(:,:,3),N)
       write(100,*)
       write(100,*)
    end if
    !-----------------------Two_Electron_Int----------------------------
    if ( N > Nmem ) then
    
        call Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)  
         
        write(*,*) " Two-electron integral complete."
        write(*,*)  
        
        if( outs == "out=4" ) then
            write(100,*) " Too many Two-electron integral !"
            write(100,*)
        end if
          
    else
    
        call Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)
        write(*,*) " Two-electron integral complete."
        write(*,*)
        
        if( outs == "out=4" ) then  
            call Two_Ele_Int_Out(N,Neri,Two_Electron,EI)
        end if
        
    end if
    !------------------------Energy between Nucleus---------------------     
    
    call Nucleus_Energy(M,Rn,MAtom,E_Nuc)
    
    !-------------------------SCF Procedure-----------------------------
    write(*,*) "=================== SCF Procedure =====================" 
    write(100,*) "=================== SCF Procedure ====================="
    write(*,*)
    write(100,*)

    H_Core =  Kinetic + N_E_Attract
    ETemp = 1.0d7
    Overlap_Temp = Overlap   
    call Mat_Sqrt(Overlap_Temp,Sqrt_S,N)
    
    !-------------------------------------------------------------------    
    if (guess == "guess=core" ) then
         call Main_Calculate_S(Sqrt_S,H_Core,MoCuA,EnA,N)
         call Main_Calculate_S(Sqrt_S,H_Core,MoCuB,EnB,N) !Initial Molecular Orbit.
    end if
    if (guess == "guess=atden" ) then
         call AtDen(M,N,Nmem,Neri,MAtom,NAtom,Overlap,Kinetic,H_core,Rn,Rt,a_Basis,c_Basis, &
                    l_Basis,bas_num,Two_Electron,EI,bas_set,t2,d2,MoCuA)        !Initial Molecular Orbit.
         MoCuB = MoCuA
    end if
    write(*,*) "                  INITIAL GUESS DONE"
    write(*,*) 
    !-------------------------------------------------------------------
    
    DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
    DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB)))

    do ii = 1,Mcyc               !The maximum number of cycles
       write(*,"(A7,I4)") " Cycle",ii
       write(100,"(A7,I4)") " Cycle",ii
       
       DensityA_T = DensityA ; DensityB_T = DensityB
       DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
       DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB))) 
       
       Temp1 = 0.0d0
       do i = 1,N 
          Temp2 = 0.0d0
          do j = 1,N
             Temp2 = Temp2 + (DensityA_T(i,j)-DensityA(i,j))**2
          end do
          Temp1 = Temp1 + Temp2
       end do
       
       Delta_DenA = dsqrt(Temp1)/real(N,8)
       
       Temp1 = 0.0d0
       do i = 1,N 
          Temp2 = 0.0d0
          do j = 1,N
             Temp2 = Temp2 + (DensityB_T(i,j)-DensityB(i,j))**2
          end do
          Temp1 = Temp1 + Temp2
       end do
       
       Delta_DenB = dsqrt(Temp1)/real(N,8)
       
       if(jud_diis) then
           if( ii < 5 ) then
               DensityA = 0.2d0*DensityA_T + 0.8d0*DensityA
               DensityB = 0.2d0*DensityB_T + 0.8d0*DensityB
           else
               DensityA = 0.5d0*DensityA_T + 0.5d0*DensityA
               DensityB = 0.5d0*DensityB_T + 0.5d0*DensityB
           end if               
       end if
       
       !------------------------ G matrix ------------------------------
       call UHF_Const_G(N,Nmem,Neri,Two_Electron,EI,DensityA,DensityB,Coul,ExcA,ExcB,t2,d2,"Y")
       !----------------------------------------------------------------
       
       idi = mod(ii,9)
       if( mod(ii,9) == 0 ) idi = 9
                  
       FockA(:,:,idi)  = H_core + Coul + ExcA
       FockB(:,:,idi)  = H_core + Coul + ExcB         !Fock Matrix
       
       if ( diis == "diis=on" ) then 
       !---------------------------  DIIS  -----------------------------
       call Error_Mat(N,FockA(:,:,idi),Overlap,DensityA,ErA(:,:,idi),emaxA)
       call Error_Mat(N,FockB(:,:,idi),Overlap,DensityB,ErB(:,:,idi),emaxB)
       
       if( ii < 5 .and. jud_diis ) ethre = 1.0d-1
       if( ii >=5 .and. jud_diis ) ethre = 5.0d-1
       
       if( max(emaxA,emaxB) < ethre) then
       
           jud_diis = .false.
           if( ii <= 9 ) then
                call UC_Coff(idi,N,ErA(:,:,1:idi),ErB(:,:,1:idi),CCoff(1:idi+1))
                FockA(:,:,0) = 0.0d0
                FockB(:,:,0) = 0.0d0
                do i = 1,idi
                    FockA(:,:,0) = FockA(:,:,0) + FockA(:,:,i)*CCoff(i)
                    FockB(:,:,0) = FockB(:,:,0) + FockB(:,:,i)*CCoff(i)
                end do 
           else
                call UC_Coff(9,N,ErA(:,:,1:9),ErB(:,:,1:9),CCoff(1:10))
                FockA(:,:,0) = 0.0d0
                FockB(:,:,0) = 0.0d0
                do i = 1,9
                    FockA(:,:,0) = FockA(:,:,0) + FockA(:,:,i)*CCoff(i)
                    FockB(:,:,0) = FockB(:,:,0) + FockB(:,:,i)*CCoff(i)
                end do 
           end if
           
       else
           FockA(:,:,0) =  FockA(:,:,idi)  
           FockB(:,:,0) =  FockB(:,:,idi)  
       end if
       !----------------------------------------------------------------  
       else
           emaxA = 0.0d0 ; emaxB = 0.0d0
           FockA(:,:,0) =  FockA(:,:,idi)  
           FockB(:,:,0) =  FockB(:,:,idi)  
       end if
                   
       call Main_Calculate_S(Sqrt_S,FockA(:,:,0),MoCuA,EnA,N)
       call Main_Calculate_S(Sqrt_S,FockB(:,:,0),MoCuB,EnB,N) !Solve FC=SCE    
        
       Temp1 = 0.0d0
       do i = 1,N
           Temp2 = 0.0d0
           do j = 1,N
               Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*H_core(j,i)
           end do
           Temp1 = Temp1 + Temp2 
       end do
       
       EeA = 0.0d0
       do i = 1,NoccA
          EeA = EeA + EnA(i)
       end do
       
       EeB = 0.0d0
       do i = 1,NoccB
          EeB = EeB + EnB(i)
       end do
       
       Ee = 0.5d0*( Temp1 + EeA + EeB )               !Calculate Energy E
       Delta_E = dabs( Ee - ETemp )
       ETemp = Ee
     
       write(*,"(A14,ES16.8)") "E(HF) = ",Ee + E_Nuc 
       write(*,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(*,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(100,"(A14,ES16.8)") "E(HF) = ",Ee + E_Nuc 
       write(100,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(100,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(*,*) 
       write(100,*)
       if( Delta_E < 10.0d0**(-Conv) .and. Delta_DenA < 1.0d-4 .and. Delta_DenB < 1.0d-4) exit
       !Judge whether the SCF procedure is converged.
       
    end do
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do i = 1,N
        Temp2 = 0.0d0
        do j = 1,N
            Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*Kinetic(j,i)
        end do
        Temp1 = Temp1 + Temp2 
    end do  
    E_T = Temp1 ; E_V = Ee + E_Nuc - E_T
    !-------------------------------------------------------------------      
    write(*,*)"======================================================="
    write(100,*)"======================================================="

    if( ii <= Mcyc-1 ) then
        write(*,"(A18,I4,A9)") " SCF Done , After",ii," cycles."
        write(100,"(A18,I4,A9)") " SCF Done , After",ii," cycles."
        write(100,"(A9,F7.4,A9,ES14.6)") "  -V/T = ",-E_V/E_T,"    KE =",E_T
        call U_Spin_S2(N,NoccA,NoccB,Overlap,MoCuA,MoCuB,S2)
        write(100,"(3(A7,F7.4))")  "<Sx>=",0.0d0,"  <Sy>=",0.0d0,"  <Sz>=",0.5d0*(NoccA-NoccB)
        write(100,"(A9,F7.4,A4,F7.4)")     "  <S**2>=",S2,"  S=",-0.5+dsqrt(1.0d0+4.0d0*S2)/2.0d0
        write(*,*)"======================================================="
       write(100,*)"======================================================="
    end if 
    
    if(  ii > Mcyc-1  ) then
        write(100,*)
        write(100,*) " Fail to convergence!"
        write(*,*)
        write(*,*) " Fail to convergence!"
        stop
    end if
    
    DensityA = matmul(MoCuA(1:N,1:NoccA),transpose(MoCuA(1:N,1:NoccA)))
    DensityB = matmul(MoCuB(1:N,1:NoccB),transpose(MoCuB(1:N,1:NoccB))) 
           
    do i = 1,N
       do j = i,N
          Temp1 = DensityA(j,i)*Overlap(j,i)
          MullikenA(i,j) = Temp1
          MullikenA(j,i) = Temp1
       end do
    end do  
    
    do i = 1,N
       do j = i,N
          Temp1 = DensityB(j,i)*Overlap(j,i)
          MullikenB(i,j) = Temp1
          MullikenB(j,i) = Temp1
       end do
    end do                   !Full Mulliken population analysis
    

    !----------------------------Print----------------------------------
    
    if( outs == "out=1" .or. outs == "out=2" .or. outs == "out=3" .or. outs == "out=4") then
       write(*,*)
       write(100,*)
       write(100,*) "Alpha Molecular Orbital Coefficients:"
       write(100,*)    
       call Mat_Out2(MoCuA,EnA,SnA,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
       write(100,*)
       write(100,*) "Beta Molecular Orbital Coefficients:"
       write(100,*)    
       call Mat_Out2(MoCuB,EnB,SnB,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
    end if
    
    if( outs == "out=2" .or. outs == "out=3" .or. outs == "out=4") then 

       !------------------------- Dipole Moment ------------------------  
       dipex = sum(DensityA*Dipole(:,:,1)) + sum(DensityB*Dipole(:,:,1))
       dipey = sum(DensityA*Dipole(:,:,2)) + sum(DensityB*Dipole(:,:,2))
       dipez = sum(DensityA*Dipole(:,:,3)) + sum(DensityB*Dipole(:,:,3))
       
       dipnx = 0.0d0 ; dipny = 0.0d0 ; dipnz = 0.0d0
       do i = 1,M
           dipnx = dipnx + MAtom(i)*Rn(1,i)
           dipny = dipny + MAtom(i)*Rn(2,i)
           dipnz = dipnz + MAtom(i)*Rn(3,i)
       end do
       
       dipx = (- dipex + dipnx)*2.54177d0
       dipy = (- dipey + dipny)*2.54177d0
       dipz = (- dipez + dipnz)*2.54177d0
       diptot = dsqrt( dipx**2 + dipy**2 + dipz**2 )
       !---------------------------------------------------------------- 
                      
       write(100,*) "Alpha Density Matrix:"
       call Mat_Out1(DensityA,N)
       write(100,*)
       write(100,*) "Beta Density Matrix:"
       call Mat_Out1(DensityB,N)
       write(100,*)
       
       write(100,*) "Full Mulliken population analysis:"
       call Mat_Out1(MullikenA+MullikenB,N)
       write(100,*)
       
       write(100,*) "Gross orbital populations:"
       write(100,*)
       call GOP_Out2(MullikenA,MullikenB,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
       
       write(100,*) "Dipole moment (field-independent basis, Debye):"
       write(100,*)
       write(100,"(3(A6,F10.4),A8,F10.4)") "X=",dipx,"Y=",dipy,"Z=",dipz,"Tot=",diptot
       write(100,*)
       
    end if
    
end subroutine UHF
!=======================================================================
!=======================================================================
subroutine ROHF(Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                NAtom,Two_Electron,EI,H_Core,MoCu,Dipole,&
                En,Ee,DensityA,DensityB,c_Basis,a_Basis,&
                l_Basis,E_Nuc,outs,bas_set,bas_num,Atom)
!------------------Open Shell，Multiplicity >= 2------------------------

    implicit none
    integer Mcyc,Conv,idi
    integer N,Nmem,Neri,NoccA,NoccB,M,bas_num(N),l_Basis(3,N)
    integer i,j,p,q,ii,MAtom(M),NAtom(N),EI(N,N)
    real(kind=8) Ee,E_Nuc,E_T,E_V,Delta_E,S2
    real(kind=8) Rn(3,M),Rt(3,N),Coul(N,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) DensityA(N,N),MullikenA(N,N),MullikenB(N,N)
    real(kind=8) DensityB(N,N),MoCu(N,N),En(N),Sqrt_S(N,N)
    real(kind=8) DensityA_T(N,N),DensityB_T(N,N)
    real(kind=8) Dipole(N,N,3),Overlap(N,N),kinetic(N,N),Fock_TT(N,N)
    real(kind=8) N_E_Attract(N,N),H_Core(N,N),Fock(N,N,0:9),Fock_T(N,N)
    real(kind=8) ErA(N,N,9),emaxA,ErB(N,N,9),emaxB,CCoff(10),ethre
    real(kind=8) Overlap_Temp(N,N),Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) dipex,dipey,dipez,dipnx,dipny,dipnz,dipx,dipy,dipz,diptot
    real(kind=8) FockA(N,N),ExcA(N,N),FockB(N,N),ExcB(N,N)
    real(kind=8) Temp1,Temp2,ETemp,Delta_DenA,Delta_DenB
    character(len=1) SnA(N),SnB(N)
    character(len=20) outs,bas_set
    character(len=20) diis,guess
    character(len=2) Atom(M),SAtom(N)
    character(len=4) Ind_Basis(N)
    character(len=4) IAtom(N)
    character(len=10) t2
    character(len=8) d2
    logical jud_diis
   
    SnA(1:NoccA) = "O" ; SnA(NoccA+1:N) = "V" 
    SnB(1:NoccB) = "O" ; SnB(NoccB+1:N) = "V" 
    jud_diis = .true.
    
    call Basis_Sets(N,M,Rn,Rt,MAtom,NAtom,bas_set,c_Basis,a_Basis,l_Basis, &
                    bas_num,Atom,IAtom,SAtom,Ind_Basis)
    
    !-----------------------One_Electron_Int----------------------------
    call Overlap_Kinetic_Int(Overlap,Kinetic,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)
    call N_E_Attractive_Int(N_E_Attract,M,N,a_Basis,c_Basis,l_Basis,bas_num,Rn,Rt,MAtom,NAtom)
    call Dipole_Int(Dipole,Overlap,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)   
     
    write(*,*) " One-electron integral complete."
     
    if( outs == "out=3" .or. outs == "out=4" ) then  
       write(100,*) "===== Overlap Matrix: ====="
       call Mat_Out1(Overlap,N)
       write(100,*)
       write(100,*)
    
       write(100,*) "===== Kinetic Matrix: ====="
       call Mat_Out1(Kinetic,N)
       write(100,*)
       write(100,*)

       write(100,*) "===== Nuclear Attraction: ====="
       call Mat_Out1(N_E_Attract,N)
       write(100,*)
       write(100,*)
       
       write(100,*) "===== Multipole matrices(X): ====="
       call Mat_Out1(Dipole(:,:,1),N)
       write(100,*)
       write(100,*)
       
       write(100,*) "===== Multipole matrices(Y): ====="
       call Mat_Out1(Dipole(:,:,2),N)
       write(100,*)
       write(100,*)
             
       write(100,*) "===== Multipole matrices(Z): ====="
       call Mat_Out1(Dipole(:,:,3),N)
       write(100,*)
       write(100,*)
    end if
    !-----------------------Two_Electron_Int----------------------------
    if ( N > Nmem ) then
    
        call Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)
         
        write(*,*) " Two-electron integral complete."
        write(*,*)  
        
        if( outs == "out=4" ) then
            write(100,*) " Too many Two-electron integral !"
            write(100,*)
        end if
          
    else
    
        call Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)
        write(*,*) " Two-electron integral complete."
        write(*,*)
        
        if( outs == "out=4" ) then  
            call Two_Ele_Int_Out(N,Neri,Two_Electron,EI)
        end if
        
    end if
    !------------------------Energy between Nucleus---------------------     

    call Nucleus_Energy(M,Rn,MAtom,E_Nuc)
    
    !-------------------------SCF Procedure-----------------------------
    write(*,*) "=================== SCF Procedure =====================" 
    write(100,*) "=================== SCF Procedure ====================="
    write(*,*)
    write(100,*)

    H_Core =  Kinetic + N_E_Attract
    Overlap_Temp = Overlap
    ETemp = 1.0d7
    call Mat_Sqrt(Overlap_Temp,Sqrt_S,N)
    
    !-------------------------------------------------------------------    
    if (guess == "guess=core" ) then
         call Main_Calculate_S(Sqrt_S,H_Core,MoCu,En,N) !Initial Molecular Orbit.
    end if
    if (guess == "guess=atden" ) then
         call AtDen(M,N,Nmem,Neri,MAtom,NAtom,Overlap,Kinetic,H_core,Rn,Rt,a_Basis,c_Basis, &
                    l_Basis,bas_num,Two_Electron,EI,bas_set,t2,d2,MoCu)        !Initial Molecular Orbit.
    end if
    write(*,*) "                  INITIAL GUESS DONE"
    write(*,*) 
    !-------------------------------------------------------------------    
    
    DensityA = matmul(MoCu(1:N,1:NoccA),transpose(MoCu(1:N,1:NoccA)))
    DensityB = matmul(MoCu(1:N,1:NoccB),transpose(MoCu(1:N,1:NoccB)))

    do ii = 1,Mcyc               !The maximum number of cycles
       write(*,"(A7,I4)") " Cycle",ii
       write(100,"(A7,I4)") " Cycle",ii

       DensityA_T = DensityA ; DensityB_T = DensityB
       DensityA = matmul(MoCu(1:N,1:NoccA),transpose(MoCu(1:N,1:NoccA)))
       DensityB = matmul(MoCu(1:N,1:NoccB),transpose(MoCu(1:N,1:NoccB)))

       Temp1 = 0.0d0
       do i = 1,N 
          Temp2 = 0.0d0
          do j = 1,N
             Temp2 = Temp2 + (DensityA_T(i,j)-DensityA(i,j))**2
          end do
          Temp1 = Temp1 + Temp2
       end do
       
       Delta_DenA = dsqrt(Temp1)/real(N,8)
       
       Temp1 = 0.0d0
       do i = 1,N 
          Temp2 = 0.0d0
          do j = 1,N
             Temp2 = Temp2 + (DensityB_T(i,j)-DensityB(i,j))**2
          end do
          Temp1 = Temp1 + Temp2
       end do
       
       Delta_DenB = dsqrt(Temp1)/real(N,8)
        
       if(jud_diis) then
           if( ii < 5 ) then
               DensityA = 0.2d0*DensityA_T + 0.8d0*DensityA
               DensityB = 0.2d0*DensityB_T + 0.8d0*DensityB
           else
               DensityA = 0.5d0*DensityA_T + 0.5d0*DensityA
               DensityB = 0.5d0*DensityB_T + 0.5d0*DensityB
           end if               
       end if
       
       !------------------------ G matrix ------------------------------
       call UHF_Const_G(N,Nmem,Neri,Two_Electron,EI,DensityA,DensityB,Coul,ExcA,ExcB,t2,d2,"Y")
       !----------------------------------------------------------------
           
       FockA = H_core + Coul + ExcA
       FockB = H_core + Coul + ExcB         !Fock Matrix
         
       do i = 1,N
         do j = i,N
           Temp1 = 0.0d0
           do p = 1,N
              Temp2 = 0.0d0
              do q = 1,N
                 Temp2 = Temp2 + MoCu(p,i)*MoCu(q,j)*FockA(p,q)
              end do
              Temp1 = Temp1 + Temp2
           end do
           Fock_T(j,i) = Temp1  
           Fock_T(i,j) = Temp1
          end do
       end do 
       
       do i = 1,N
         do j = i,N
           Temp1 = 0.0d0
           do p = 1,N
              Temp2 = 0.0d0
              do q = 1,N
                 Temp2 = Temp2 + MoCu(p,i)*MoCu(q,j)*FockB(p,q)
              end do
              Temp1 = Temp1 + Temp2
           end do
           Fock_TT(j,i) = Temp1  
           Fock_TT(i,j) = Temp1
          end do
       end do 
       
       idi = mod(ii,9)
       if( mod(ii,9) == 0 ) idi = 9   
             
       call ROHF_Ham(Fock(:,:,idi),Fock_T,Fock_TT,SnA,SnB,N)
       
       call Inverse(transpose(MoCu),Fock_T,N)
       call Inverse(MoCu,Fock_TT,N)
       Fock(:,:,idi) = matmul(matmul(Fock_T,Fock(:,:,idi)),Fock_TT)
       
       if ( diis == "diis=on" ) then       
       !---------------------------  DIIS  -----------------------------
       call Error_Mat(N,Fock(:,:,idi),Overlap,DensityA,ErA(:,:,idi),emaxA)
       call Error_Mat(N,Fock(:,:,idi),Overlap,DensityB,ErB(:,:,idi),emaxB)
       
       if( ii < 5 .and. jud_diis ) ethre = 1.0d-1
       if( ii >=5 .and. jud_diis ) ethre = 5.0d-1
       
       if( max(emaxA,emaxB) < ethre ) then
       
           jud_diis = .false.
           if( ii <= 9 ) then
                call UC_Coff(idi,N,ErA(:,:,1:idi),ErB(:,:,1:idi),CCoff(1:idi+1))
                Fock(:,:,0) = 0.0d0
                do i = 1,idi
                   Fock(:,:,0) = Fock(:,:,0) + Fock(:,:,i)*CCoff(i)
                end do 
           else
                call UC_Coff(9,N,ErA(:,:,1:9),ErB(:,:,1:9),CCoff(1:10))
                Fock(:,:,0) = 0.0d0
                do i = 1,9
                   Fock(:,:,0) = Fock(:,:,0) + Fock(:,:,i)*CCoff(i)
                end do 
           end if
           
       else
           Fock(:,:,0) =  Fock(:,:,idi)    
       end if
       !----------------------------------------------------------------
       else
           emaxA = 0.0d0 ; emaxB = 0.0d0
           Fock(:,:,0) =  Fock(:,:,idi)    
       end if
                 
       call Main_Calculate_S(Sqrt_S,Fock(:,:,0),MoCu,En,N)      
        
       Temp1 = 0.0d0
       do i = 1,N
           Temp2 = 0.0d0
           do j = 1,N
               Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*H_core(j,i) + &
                       DensityA(j,i)*FockA(j,i) + DensityB(j,i)*FockB(j,i) 
           end do
           Temp1 = Temp1 + Temp2 
       end do

       Ee = Temp1*0.5d0                !Calculate Energy E
       Delta_E = dabs( Ee - ETemp )
       
       ETemp = Ee
       write(*,"(A14,ES16.8)") "E(HF) = ",Ee + E_Nuc 
       write(*,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(*,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(100,"(A14,ES16.8)") "E(HF) = ",Ee + E_Nuc 
       write(100,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(100,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(*,*) 
       write(100,*)
       if( Delta_E < 10.0d0**(-Conv) .and. Delta_DenA < 1.0d-4 .and. Delta_DenB < 1.0d-4) exit
       !Judge whether the SCF procedure is converged.
       
    end do
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do i = 1,N
        Temp2 = 0.0d0
        do j = 1,N
            Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*Kinetic(j,i)
        end do
        Temp1 = Temp1 + Temp2 
    end do  
    E_T = Temp1 ; E_V = Ee + E_Nuc - E_T
    !-------------------------------------------------------------------      
    write(*,*)"======================================================="
    write(100,*)"=======================================================" 
    
    if( ii <= Mcyc-1 ) then
        write(*,"(A18,I4,A9)") " SCF Done , After",ii," cycles."
        write(100,"(A18,I4,A9)") " SCF Done , After",ii," cycles."
        write(100,"(A9,F7.4,A9,ES14.6)") "  -V/T = ",-E_V/E_T,"    KE =",E_T
        S2= 0.5d0*(NoccA-NoccB)*(0.5d0*(NoccA-NoccB)+1)
        write(100,"(3(A7,F7.4))")  "<Sx>=",0.0d0,"  <Sy>=",0.0d0,"  <Sz>=",0.5d0*(NoccA-NoccB)
        write(100,"(A9,F7.4,A4,F7.4)")     "  <S**2>=",S2,"  S=",-0.5+dsqrt(1.0d0+4.0d0*S2)/2.0d0
        write(*,*)"======================================================="
       write(100,*)"======================================================="
    end if 
    
    if(  ii > Mcyc-1  ) then
        write(100,*)
        write(100,*) " Fail to convergence!"
        write(*,*)
        write(*,*) " Fail to convergence!"
        stop
    end if
    
    DensityA = matmul(MoCu(1:N,1:NoccA),transpose(MoCu(1:N,1:NoccA)))
    DensityB = matmul(MoCu(1:N,1:NoccB),transpose(MoCu(1:N,1:NoccB)))  
           
    do i = 1,N
       do j = i,N
          Temp1 = DensityA(j,i)*Overlap(j,i)
          MullikenA(i,j) = Temp1
          MullikenA(j,i) = Temp1
       end do
    end do  
    
    do i = 1,N
       do j = i,N
          Temp1 = DensityB(j,i)*Overlap(j,i)
          MullikenB(i,j) = Temp1
          MullikenB(j,i) = Temp1
       end do
    end do                     !Full Mulliken population analysis
    
    !----------------------------Print----------------------------------
    
    if( outs == "out=1" .or. outs == "out=2" .or. outs == "out=3" .or. outs == "out=4") then
       write(*,*)
       write(100,*)
       write(100,*) "Molecular Orbital Coefficients:"
       write(100,*)    
       call Mat_Out2(MoCu,En,SnA,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
    end if
    
    if( outs == "out=2" .or. outs == "out=3" .or. outs == "out=4") then 
    
       !------------------------- Dipole Moment ------------------------  
       dipex = sum(DensityA*Dipole(:,:,1)) + sum(DensityB*Dipole(:,:,1))
       dipey = sum(DensityA*Dipole(:,:,2)) + sum(DensityB*Dipole(:,:,2))
       dipez = sum(DensityA*Dipole(:,:,3)) + sum(DensityB*Dipole(:,:,3))
       
       dipnx = 0.0d0 ; dipny = 0.0d0 ; dipnz = 0.0d0
       do i = 1,M
           dipnx = dipnx + MAtom(i)*Rn(1,i)
           dipny = dipny + MAtom(i)*Rn(2,i)
           dipnz = dipnz + MAtom(i)*Rn(3,i)
       end do
       
       dipx = (- dipex + dipnx)*2.54177d0
       dipy = (- dipey + dipny)*2.54177d0
       dipz = (- dipez + dipnz)*2.54177d0
       diptot = dsqrt( dipx**2 + dipy**2 + dipz**2 )
       !----------------------------------------------------------------       
                  
       write(100,*) "Alpha Density Matrix:"
       call Mat_Out1(DensityA,N)
       write(100,*)
       write(100,*) "Beta Density Matrix:"
       call Mat_Out1(DensityB,N)
       write(100,*)
       
       write(100,*) "Full Mulliken population analysis:"
       call Mat_Out1(MullikenA+MullikenB,N)
       write(100,*)
       
       write(100,*) "Gross orbital populations:"
       write(100,*)
       call GOP_Out2(MullikenA,MullikenB,N,IAtom,SAtom,Ind_Basis)
       write(100,*)
       
       write(100,*) "Dipole moment (field-independent basis, Debye):"
       write(100,*)
       write(100,"(3(A6,F10.4),A8,F10.4)") "X=",dipx,"Y=",dipy,"Z=",dipz,"Tot=",diptot
       write(100,*)
       
    end if
    
end subroutine ROHF
!=======================================================================
subroutine ROHF_Ham(Fock,FockA,FockB,SnA,SnB,N)

    implicit None
    integer  N,i,j
    real(kind=8) Fock(N,N),FockA(N,N),FockB(N,N)
    real(kind=8) ac,bc,ao,bo,av,bv,fij,pt
    character(len=1) SnA(N),SnB(N)
    logical  ic,io,iv,jc,jo,jv

    ac = 1.0d0/3.0d0 ; bc = 2.0d0/3.0d0
    ao = 1.0d0/3.0d0 ; bo = 1.0d0/3.0d0
    av = 2.0d0/3.0d0 ; bv = 1.0d0/3.0d0
    pt = 0.5d0
    
    do j = 1,N
       jc = SnA(j) == "O" .and. SnB(j) == "O" 
       jo = SnA(j) == "O" .and. SnB(j) /= "O" 
       jv = SnA(j) /= "O" .and. SnB(j) /= "O" 
       do i = 1,j
          ic = SnA(i) == "O" .and. SnB(i) == "O" 
          io = SnA(i) == "O" .and. SnB(i) /= "O" 
          iv = SnA(i) /= "O" .and. SnB(i) /= "O" 

          if (ic .and. jc) fij = ac*FockA(i,j) + bc*FockB(i,j)
          if (ic .and. jo) fij = FockB(i,j)
          if (ic .and. jv) fij = pt*(FockA(i,j)+FockB(i,j))
          if (io .and. jo) fij = ao*FockA(i,j) + bo*FockB(i,j)
          if (io .and. jv) fij = FockA(i,j)
          if (iv .and. jv) fij = av*FockA(i,j) + bv*FockB(i,j)
          Fock(i, j) = fij
          Fock(j, i) = fij
       end do
  end do
  
end subroutine ROHF_Ham
!=======================================================================


