subroutine DFT_Close(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,Nocc,M,Rn,Rt,MAtom,NAtom, &
                     Two_Electron,EI,H_Core,MoCu,En,Ee,Density,Dipole,&
                     c_Basis,a_Basis,l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
!------------------Close Shell，Multiplicity = 1------------------------
!-----------------------------------------------------------------------
!                     E = T + E_Nuc + Coul + Exc                       !
!    Vxc = delta{Exc}/delta{rho} --> Vrho = delta{Func_E}/delta{rho}   !
!                 Func_E --> Exc  ;  Vrho --> Vxc ;                    !
!-----------------------------------------------------------------------

    implicit none
    integer Node,Node1
    integer idi,Mcyc,Conv
    integer N,Nmem,Neri,Nocc,M,bas_num(N),MAtom(M),NAtom(N)
    integer La(3),EI(N,N)
    integer i,j,k,l,ii,i1,k1,end_pot(M),l_Basis(3,N)
    real(kind=8) Ee,E_Nuc,E_T,E_V,E_Exc,Delta_E,E_XC
    real(kind=8) Rn(3,M),Rt(3,N),a,Ra(3)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) MoCu(N,N),En(N),Density(N,N),Mulliken(N,N)
    real(kind=8) Overlap(N,N),kinetic(N,N),Sqrt_S(N,N)
    real(kind=8) N_E_Attract(N,N),H_Core(N,N),Dipole(N,N,3)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2),Density_T(N,N)
    real(kind=8) Fock(N,N,0:7),Coul(N,N),Exc(N,N),XC(N,N)
    real(kind=8) Er(N,N,7),emax,CCoff(8),ethre
    real(kind=8) gridx(Node,M),gridy(Node,M),gridz(Node,M)
    real(kind=8) W_i(Node,M),Func_E(Node,M),Vrho(Node,M),Vsigma(Node,M)
    real(kind=8) dipex,dipey,dipez,dipnx,dipny,dipnz,dipx,dipy,dipz,diptot
    real(kind=8) Temp1,Temp2,ETemp,Overlap_Temp(N,N),w0,Delta_Den,C0
    real(kind=8) D_Phi(3,N),ATemp0,Phi(N),x,y,z
    real(kind=8) A3Temp0(3),Rho_grid(Node,M,3)
    character(len=1) Sn(N)
    character(len=20) outs,bas_set,method
    character(len=20) diis,guess
    character(len=4) IAtom(N)
    character(len=4) Ind_Basis(N)
    character(len=2) Atom(M),SAtom(N)
    character(len=10) t2
    character(len=8) d2
    character(len=1) CE
    logical jud_diis,jud_node
    
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
                    l_Basis,bas_num,Two_Electron,EI,bas_set,t2,d2,MoCu)        !Initial Molecular Orbit.
    end if
    write(*,*) "                  INITIAL GUESS DONE"
    write(*,*) 
    !-------------------------------------------------------------------
    Density = matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
    
    if(Node == 22650) Node1 = 8760
    if(Node == 8760) Node1 = 2960
    jud_node = .true.
    call Int_Dep(Node1,M,Rn,MAtom,gridx(1:Node1,1:M),gridy(1:Node1,1:M),gridz(1:Node1,1:M), &
                 W_i(1:Node1,1:M),end_pot)
    !-------------------------------------------------------------------
    if( method == "b3lyp" .or. method == "B3LYP" .or. method == "pbe1pbe" &
          .or. method == "PBE1PBE" .or. method == "b3pw91" .or. method == "B3PW91") then
          C0 = 0.3d0
          CE = "Y"
    else
          C0 = 0.5d0
          CE = "C"
    end if
    !-------------------------------------------------------------------
    do ii = 1,Mcyc                 !The maximum number of cycles
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
       
       if(jud_diis) Density = C0*Density_T + (1.0d0-C0)*Density
       
       !------------------------ G matrix ------------------------------
       call RDFT_Const_G(N,Nmem,Neri,Two_Electron,EI,Density,Coul,Exc,t2,d2,CE)
       !----------------------------------------------------------------
       
       !------------------- rho,Func_E,Vrho,Vsigma ---------------------
       call rho_Func_Close(Node1,M,Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,method, &
                   gridx(1:Node1,1:M),gridy(1:Node1,1:M),gridz(1:Node1,1:M),end_pot,Rho_grid(1:Node1,1:M,1:3), &
                   Func_E(1:Node1,1:M),Vrho(1:Node1,1:M),Vsigma(1:Node1,1:M))   

       !------------------------- XC Matrix ----------------------------
       XC = 0.0d0
       if( method == "XAlpha" .or. method == "xalpha" &
           .or. method == "LSDA" .or. method == "lsda" ) then

         do i = 1,M
            do j = 1,end_pot(i)    
              x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)  
              !---------------------------------------------------------
              do i1 = 1,N              
                 Phi(i1) = 0.0d0
                 D_Phi(:,i1) = 0.0d0
                 do k1 = 1,bas_num(i1)
                      a = a_Basis(k1,i1,NAtom(i1)) ; Ra = Rt(:,i1) 
                      La = l_Basis(:,i1)
                      call Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,"N")
                      Phi(i1) = Phi(i1) + ATemp0*c_Basis(k1,i1,NAtom(i1))
                 end do
               end do             
               !--------------------------------------------------------
               do k = 1,N
                  do l = k,N
                  
                      XC(l,k) =  XC(l,k) + Phi(k)*Phi(l)*Vrho(j,i)*W_i(j,i)

                  end do       
               end do 
               !--------------------------------------------------------
            end do
         end do   
          
       else
      
         do i = 1,M
            do j = 1,end_pot(i)  
              x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)  
              !---------------------------------------------------------
              do i1 = 1,N              
                 Phi(i1) = 0.0d0
                 D_Phi(:,i1) = 0.0d0
                 do k1 = 1,bas_num(i1)
                      a = a_Basis(k1,i1,NAtom(i1)) ; Ra = Rt(:,i1) 
                      La = l_Basis(:,i1)
                      call Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,"Y")
                      Phi(i1) = Phi(i1) + ATemp0*c_Basis(k1,i1,NAtom(i1))
                      D_Phi(:,i1) = D_Phi(:,i1) + A3Temp0*c_Basis(k1,i1,NAtom(i1))
                 end do
               end do             
               !---------------------------------------------------------    
               do k = 1,N
                  do l = k,N
                      
                        Temp1 = Phi(k)*Phi(l)*Vrho(j,i)
                        Temp2 = 2.0d0*Vsigma(j,i)*sum( Rho_grid(j,i,:)*( Phi(k)*D_Phi(:,l) + Phi(l)*D_Phi(:,k) ) )
                        XC(l,k) =  XC(l,k) + ( Temp1 + Temp2 )*W_i(j,i)
                        
                  end do       
               end do 
               !--------------------------------------------------------
            end do
         end do
        
       end if 
        
       do k = 1,N
          do l = k,N
              XC(k,l) = XC(l,k)  
          end do
       end do                
       !----------------------------------------------------------------
       E_Exc = 0.0d0 ; w0 = 0.0d0
       if( method == "b3lyp" .or. method == "B3LYP" .or. method == "pbe1pbe" &
          .or. method == "PBE1PBE" .or. method == "b3pw91" .or. method == "B3PW91") then
          
          if( method == "b3pw91" .or. method == "B3PW91" ) w0 = 0.2d0
          if( method == "b3lyp" .or. method == "B3LYP" ) w0 = 0.2d0
          if( method == "pbe1pbe" .or. method == "PBE1PBE" ) w0 = 0.25d0      
         
          XC = XC + w0*Exc 
          
          Temp1 = 0.0d0
          do i = 1,N
              Temp2 = 0.0d0
              do j = 1,N
                   Temp2 = Temp2 + Density(j,i)*Exc(j,i)
              end do
              Temp1 = Temp1 + Temp2 
          end do     
          E_Exc = Temp1
               
       end if
       !----------------------------------------------------------------
       
       idi = mod(ii,7)
       if( mod(ii,7) == 0 ) idi = 7   
           
       Fock(:,:,idi) = H_core + Coul + XC     !Fock Matrix
       
       if ( diis == "diis=on" ) then        
       !---------------------------  DIIS  -----------------------------
       call Error_Mat(N,Fock(:,:,idi),Overlap,Density,Er(:,:,idi),emax)
       
       ethre = 5.0d-1
       
       if( emax < ethre) then
       
           jud_diis = .false.
           if( ii <= 7 ) then
                call C_Coff(idi,N,Er(:,:,1:idi),CCoff(1:idi+1))
                Fock(:,:,0) = 0.0d0
                do i = 1,idi
                   Fock(:,:,0) = Fock(:,:,0) + Fock(:,:,i)*CCoff(i)
                end do 
           else
                call C_Coff(7,N,Er(:,:,1:7),CCoff(1:8))
                Fock(:,:,0) = 0.0d0
                do i = 1,7
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
                Temp2 = Temp2 + Density(j,i)*(2.0d0*H_core(j,i) + Coul(j,i))
           end do
           Temp1 = Temp1 + Temp2 
       end do
       
       call Int_Func_E_XC(Node1,M,Func_E(1:Node1,1:M),W_i(1:Node1,1:M),end_pot,E_XC)
       Ee = Temp1 + E_XC + w0*E_Exc               !Calculate Energy E
       Delta_E = dabs( Ee - ETemp )
       
       write(*,"(A14,ES16.8)") "E(DFT)= ",Ee + E_Nuc 
       write(*,"(A14,ES14.6)") "DIIS_ER=",emax
       write(*,"(A13,ES15.6,A10,ES16.6)") "DDen  =",Delta_Den,"DE =",Delta_E
       write(100,"(A14,ES16.8)") "E(DFT)= ",Ee + E_Nuc 
       write(100,"(A14,ES14.6)") "DIIS_ER=",emax
       write(100,"(A13,ES15.6,A10,ES16.6)") "DDen  =",Delta_Den,"DE =",Delta_E
       write(*,*)
       write(100,*)
       if( Delta_E < 10.0d0**(-Conv) .and. Delta_Den < 1.0d-4 ) exit 
       !Judge whether the SCF procedure is converged.
       
       ETemp = Ee
       
       if( Delta_E < 10.0d0**(-5.0d0) .and. jud_node) then 
           jud_node = .false.
           Node1 = Node
           write(*,*) "     Initial convergence to 1.0D-05 achieved. "
           write(*,*) "     Increase integral accuracy."
           write(*,*)
           write(100,*) "     Initial convergence to 1.0D-05 achieved. "
           write(100,*) "     Increase integral accuracy."
           write(100,*)
           call Int_Dep(Node1,M,Rn,MAtom,gridx(1:Node1,1:M),gridy(1:Node1,1:M), &
                gridz(1:Node1,1:M), W_i(1:Node1,1:M),end_pot)
       end if
       
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
    
    if(  ii > Mcyc-1  ) then
        write(100,*)
        write(100,*) " Fail to convergence!"
        write(*,*)
        write(*,*) " Fail to convergence!"
        stop
    end if
       
    Density = matmul(MoCu(1:N,1:Nocc),transpose(MoCu(1:N,1:Nocc)))
    Density = 2.0d0*Density   !Density Matrix
    
    do i = 1,N
       do j = i,N
          Temp1 = Density(j,i)*Overlap(j,i)
          Mulliken(i,j) = Temp1
          Mulliken(j,i) = Temp1
       end do
    end do                    !Full Mulliken population analysis
    
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
    
end subroutine DFT_Close
!=======================================================================
!=======================================================================
subroutine DFT_Open(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                    NAtom,Two_Electron,EI,H_Core,MoCuA,MoCuB,Dipole,&
                    EnA,EnB,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                    l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
!------------------Open Shell，Multiplicity >= 2------------------------
!-----------------------------------------------------------------------
!                     E = T + E_Nuc + Coul + Exc                       !
!    Vxc = delta{Exc}/delta{rho} --> Vrho = delta{Func_E}/delta{rho}   !
!                 Func_E --> Exc  ;  Vrho --> Vxc ;                    !
!-----------------------------------------------------------------------

    implicit none
    integer Node,Node1
    integer Mcyc,Conv,idi
    integer N,Nmem,Neri,NoccA,NoccB,M,bas_num(N),MAtom(M),NAtom(N)
    integer La(3),EI(N,N)
    integer i,j,k,l,i1,k1,ii,end_pot(M),l_Basis(3,N)
    real(kind=8) Ee,E_Nuc,E_T,E_V,E_Exc,Delta_E,E_XC,S2
    real(kind=8) Rn(3,M),Rt(3,N),a,Ra(3)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) MoCuA(N,N),EnA(N),DensityA(N,N),MullikenA(N,N)
    real(kind=8) MoCuB(N,N),EnB(N),DensityB(N,N),MullikenB(N,N)
    real(kind=8) Overlap(N,N),kinetic(N,N),Sqrt_S(N,N),Dipole(N,N,3)
    real(kind=8) DensityA_T(N,N),DensityB_T(N,N)
    real(kind=8) N_E_Attract(N,N),H_Core(N,N),ExcA(N,N),ExcB(N,N)
    real(kind=8) Overlap_Temp(N,N),Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) FockA(N,N,0:7),FockB(N,N,0:7),Coul(N,N),XCA(N,N),XCB(N,N)
    real(kind=8) ErA(N,N,7),emaxA,CCoff(8)
    real(kind=8) ErB(N,N,7),emaxB,ethre
    real(kind=8) gridx(Node,M),gridy(Node,M),gridz(Node,M),VsigmaA(Node,M),VsigmaB(Node,M)
    real(kind=8) dipex,dipey,dipez,dipnx,dipny,dipnz,dipx,dipy,dipz,diptot
    real(kind=8) W_i(Node,M),VrhoA(Node,M),VrhoB(Node,M),Func_E(Node,M),VsigmaAB(Node,M)
    real(kind=8) A3Temp0(3),ATemp0,x,y,z
    real(kind=8) Phi(N),D_Phi(3,N)
    real(kind=8) Rho_gridA(Node,M,3),Rho_gridB(Node,M,3)
    real(kind=8) Temp1,Temp2,ETemp,w0,Delta_DenA,Delta_DenB,C0
    character(len=1) SnA(N),SnB(N)
    character(len=20) outs,bas_set,method
    character(len=20) diis,guess
    character(len=4) IAtom(N)
    character(len=2) Atom(M),SAtom(N)
    character(len=4) Ind_Basis(N)
    character(len=10) t2
    character(len=8) d2
    character(len=1) CE
    logical jud_diis,jud_node
   
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
    
    if(Node == 22650) Node1 = 8760
    if(Node == 8760) Node1 = 2960
    jud_node = .true.
    call Int_Dep(Node1,M,Rn,MAtom,gridx(1:Node1,1:M),gridy(1:Node1,1:M),gridz(1:Node1,1:M), &
                 W_i(1:Node1,1:M),end_pot)
    !-------------------------------------------------------------------
    if( method == "b3lyp" .or. method == "B3LYP" .or. method == "pbe1pbe" &
          .or. method == "PBE1PBE" .or. method == "b3pw91" .or. method == "B3PW91") then
          C0 = 0.3d0
          CE = "Y"
    else
          C0 = 0.5d0
          CE = "C"
    end if
    !-------------------------------------------------------------------
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
           DensityA = C0*DensityA_T + (1.0d0-C0)*DensityA
           DensityB = C0*DensityB_T + (1.0d0-C0)*DensityB
       end if
       
       !------------------------ G matrix ------------------------------
       call UHF_Const_G(N,Nmem,Neri,Two_Electron,EI,DensityA,DensityB,Coul,ExcA,ExcB,t2,d2,CE)
       !----------------------------------------------------------------
       
       !------------------- rho,Func_E,Vrho,Vsigma ---------------------
       call rho_Func_Open(Node1,M,DensityA,DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,method, &
                          gridx(1:Node1,1:M),gridy(1:Node1,1:M),gridz(1:Node1,1:M), &
                          end_pot,Rho_gridA(1:Node1,1:M,1:3),Rho_gridB(1:Node1,1:M,1:3), &
                          Func_E(1:Node1,1:M),VrhoA(1:Node1,1:M),VrhoB(1:Node1,1:M), &
                          VsigmaA(1:Node1,1:M),VsigmaB(1:Node1,1:M),VsigmaAB(1:Node1,1:M)) 
     
       !------------------------- XC Matrix ----------------------------       
       XCA = 0.0d0 ; XCB = 0.0d0
       if( method == "XAlpha" .or. method == "xalpha" &
           .or. method == "LSDA" .or. method == "lsda" ) then
      
         do i = 1,M
            do j = 1,end_pot(i)    
              x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)  
              !---------------------------------------------------------
              do i1 = 1,N              
                 Phi(i1) = 0.0d0
                 D_Phi(:,i1) = 0.0d0
                 do k1 = 1,bas_num(i1)
                      a = a_Basis(k1,i1,NAtom(i1)) ; Ra = Rt(:,i1) 
                      La = l_Basis(:,i1)
                      call Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,"N")
                      Phi(i1) = Phi(i1) + ATemp0*c_Basis(k1,i1,NAtom(i1))
                 end do
               end do             
               !--------------------------------------------------------
               do k = 1,N
                  do l = k,N
                  
                      XCA(l,k) =  XCA(l,k) + Phi(k)*Phi(l)*VrhoA(j,i)*W_i(j,i)
                      XCB(l,k) =  XCB(l,k) + Phi(k)*Phi(l)*VrhoB(j,i)*W_i(j,i)

                  end do       
               end do 
               !--------------------------------------------------------
            end do
         end do 
       
       else
      
         do i = 1,M
            do j = 1,end_pot(i)  
              x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)  
              !---------------------------------------------------------
              do i1 = 1,N              
                 Phi(i1) = 0.0d0
                 D_Phi(:,i1) = 0.0d0
                 do k1 = 1,bas_num(i1)
                      a = a_Basis(k1,i1,NAtom(i1)) ; Ra = Rt(:,i1) 
                      La = l_Basis(:,i1)
                      call Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,"Y")
                      Phi(i1) = Phi(i1) + ATemp0*c_Basis(k1,i1,NAtom(i1))
                      D_Phi(:,i1) = D_Phi(:,i1) + A3Temp0*c_Basis(k1,i1,NAtom(i1))
                 end do
               end do             
               !---------------------------------------------------------    
               do k = 1,N
                  do l = k,N
                      
                      Temp1 = Phi(k)*Phi(l)*VrhoA(j,i) 
                      Temp2 = sum( ( 2.0d0*VsigmaA(j,i)*Rho_gridA(j,i,:) + VsigmaAB(j,i)*Rho_gridB(j,i,:) )*&
                              ( Phi(k)*D_Phi(:,l) + Phi(l)*D_Phi(:,k) ) )
                      XCA(l,k) =  XCA(l,k) + ( Temp1 + Temp2 )*W_i(j,i)
                      
                      Temp1 = Phi(k)*Phi(l)*VrhoB(j,i) 
                      Temp2 = sum( ( 2.0d0*VsigmaB(j,i)*Rho_gridB(j,i,:) + VsigmaAB(j,i)*Rho_gridA(j,i,:) )*&
                              ( Phi(k)*D_Phi(:,l) + Phi(l)*D_Phi(:,k) ) )
                      XCB(l,k) =  XCB(l,k) + ( Temp1 + Temp2 )*W_i(j,i)
                                  
                  end do       
               end do 
               !--------------------------------------------------------
            end do
         end do  
       
       end if 
       
       do k = 1,N
          do l = k,N
              XCA(k,l) = XCA(l,k) 
              XCB(k,l) = XCB(l,k)  
          end do
       end do                    
       !----------------------------------------------------------------
       E_Exc = 0.0d0 ; w0 = 0.0d0
       if( method == "b3lyp" .or. method == "B3LYP" .or. method == "pbe1pbe" &
           .or. method == "PBE1PBE" .or. method == "b3pw91" .or. method == "B3PW91") then
           
          if( method == "b3lyp" .or. method == "B3LYP" ) w0 = 0.2d0
          if( method == "b3pw91" .or. method == "B3PW91" ) w0 = 0.2d0
          if( method == "pbe1pbe" .or. method == "PBE1PBE" ) w0 = 0.25d0
         
          XCA = XCA + w0*ExcA
          XCB = XCB + w0*ExcB 
          
          Temp1 = 0.0d0
          do i = 1,N
              Temp2 = 0.0d0
              do j = 1,N
                   Temp2 = Temp2 + DensityA(j,i)*ExcA(j,i)+DensityB(j,i)*ExcB(j,i)
              end do
              Temp1 = Temp1 + Temp2 
          end do     
          E_Exc = Temp1*0.5d0
               
       end if
       !----------------------------------------------------------------
       
       idi = mod(ii,7)
       if( mod(ii,7) == 0 ) idi = 7
                                           
       FockA(:,:,idi) = H_core + Coul + XCA
       FockB(:,:,idi) = H_core + Coul + XCB        !Fock Matrix
       
       if ( diis == "diis=on" ) then 
       !---------------------------  DIIS  -----------------------------
       call Error_Mat(N,FockA(:,:,idi),Overlap,DensityA,ErA(:,:,idi),emaxA)
       call Error_Mat(N,FockB(:,:,idi),Overlap,DensityB,ErB(:,:,idi),emaxB)
       
       ethre = 5.0d-1

       if( max(emaxA,emaxB) < ethre) then
       
           jud_diis = .false.
           if( ii <= 7 ) then
                call UC_Coff(idi,N,ErA(:,:,1:idi),ErB(:,:,1:idi),CCoff(1:idi+1))
                FockA(:,:,0) = 0.0d0
                FockB(:,:,0) = 0.0d0
                do i = 1,idi
                    FockA(:,:,0) = FockA(:,:,0) + FockA(:,:,i)*CCoff(i)
                    FockB(:,:,0) = FockB(:,:,0) + FockB(:,:,i)*CCoff(i)
                end do 
           else
                call UC_Coff(7,N,ErA(:,:,1:7),ErB(:,:,1:7),CCoff(1:8))
                FockA(:,:,0) = 0.0d0
                FockB(:,:,0) = 0.0d0
                do i = 1,7
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
                Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*(H_core(j,i) + 0.5d0*Coul(j,i))
           end do
           Temp1 = Temp1 + Temp2 
       end do
       
       call Int_Func_E_XC(Node1,M,Func_E(1:Node1,1:M),W_i(1:Node1,1:M),end_pot,E_XC)
             
       Ee = Temp1 + E_XC + w0*E_Exc           !Calculate Energy E
       Delta_E = dabs( Ee - ETemp )
       ETemp = Ee
       
       write(*,"(A14,ES16.8)") "E(DFT)= ",Ee + E_Nuc 
       write(*,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(*,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(100,"(A14,ES16.8)") "E(DFT)= ",Ee + E_Nuc 
       write(100,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(100,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(*,*) 
       write(100,*)
       if( Delta_E < 10.0d0**(-Conv) .and. Delta_DenA < 1.0d-4 .and. Delta_DenB < 1.0d-4) exit
       !Judge whether the SCF procedure is converged.
       
       if( Delta_E < 10.0d0**(-5.0d0) .and. jud_node) then 
           jud_node = .false.
           Node1 = Node
           write(*,*) "     Initial convergence to 1.0D-05 achieved. "
           write(*,*) "     Increase integral accuracy."
           write(*,*)
           write(100,*) "     Initial convergence to 1.0D-05 achieved. "
           write(100,*) "     Increase integral accuracy."
           write(100,*)
           call Int_Dep(Node1,M,Rn,MAtom,gridx(1:Node1,1:M),gridy(1:Node1,1:M), &
                gridz(1:Node1,1:M), W_i(1:Node1,1:M),end_pot)
       end if
              
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
    end do                     !Full Mulliken population analysis
    
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
    
end subroutine DFT_Open
!=======================================================================
!=======================================================================
subroutine DFT_RO(Node,Mcyc,Conv,diis,guess,N,Nmem,Neri,NoccA,NoccB,M,Rn,Rt,MAtom, &
                    NAtom,Two_Electron,EI,H_Core,MoCu,Dipole, &
                    En,Ee,DensityA,DensityB,c_Basis,a_Basis, &
                    l_Basis,E_Nuc,outs,bas_set,method,bas_num,Atom)
!------------------Open Shell，Multiplicity >= 2------------------------
!-----------------------------------------------------------------------
!                     E = T + E_Nuc + Coul + Exc                       !
!    Vxc = delta{Exc}/delta{rho} --> Vrho = delta{Func_E}/delta{rho}   !
!                 Func_E --> Exc  ;  Vrho --> Vxc ;                    !
!-----------------------------------------------------------------------

    implicit none
    integer Node,Node1
    integer Mcyc,Conv,idi
    integer N,Nmem,Neri,NoccA,NoccB,M,bas_num(N),MAtom(M),NAtom(N)
    integer La(3),EI(N,N)
    integer i,j,k,l,i1,k1,ii,p,q,end_pot(M),l_Basis(3,N)
    real(kind=8) Ee,E_Nuc,E_T,E_V,E_Exc,Delta_E,E_XC,S2
    real(kind=8) Rn(3,M),Rt(3,N),a,Ra(3)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) MoCu(N,N),En(N),DensityA(N,N),MullikenA(N,N)
    real(kind=8) DensityB(N,N),MullikenB(N,N)
    real(kind=8) Overlap(N,N),kinetic(N,N),Sqrt_S(N,N)
    real(kind=8) DensityA_T(N,N),DensityB_T(N,N),Dipole(N,N,3)
    real(kind=8) N_E_Attract(N,N),H_Core(N,N),ExcA(N,N),ExcB(N,N)
    real(kind=8) Overlap_Temp(N,N),Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) FockA(N,N),FockB(N,N),Coul(N,N),XCA(N,N),XCB(N,N)
    real(kind=8) ErA(N,N,7),emaxA,CCoff(8),Fock(N,N,0:7)
    real(kind=8) ErB(N,N,7),emaxB,ethre,Fock_T(N,N),Fock_TT(N,N)
    real(kind=8) gridx(Node,M),gridy(Node,M),gridz(Node,M),VsigmaA(Node,M),VsigmaB(Node,M)
    real(kind=8) dipex,dipey,dipez,dipnx,dipny,dipnz,dipx,dipy,dipz,diptot
    real(kind=8) W_i(Node,M),VrhoA(Node,M),VrhoB(Node,M),Func_E(Node,M),VsigmaAB(Node,M)
    real(kind=8) A3Temp0(3),ATemp0,x,y,z
    real(kind=8) Phi(N),D_Phi(3,N)
    real(kind=8) Rho_gridA(Node,M,3),Rho_gridB(Node,M,3)
    real(kind=8) Temp1,Temp2,ETemp,w0,Delta_DenA,Delta_DenB,C0
    character(len=1) SnA(N),SnB(N)
    character(len=20) outs,bas_set,method
    character(len=20) diis,guess
    character(len=4) IAtom(N)
    character(len=2) Atom(M),SAtom(N)
    character(len=4) Ind_Basis(N)
    character(len=10) t2
    character(len=8) d2
    character(len=1) CE
    logical jud_diis,jud_node
   
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
    
    if(Node == 22650) Node1 = 8760
    if(Node == 8760) Node1 = 2960
    jud_node = .true.
    call Int_Dep(Node1,M,Rn,MAtom,gridx(1:Node1,1:M),gridy(1:Node1,1:M),gridz(1:Node1,1:M), &
                 W_i(1:Node1,1:M),end_pot)
    !-------------------------------------------------------------------
    if( method == "rob3lyp" .or. method == "ropbe1pbe" .or. method == "rob3pw91") then
          C0 = 0.3d0
          CE = "Y"
    else
          C0 = 0.5d0
          CE = "C"
    end if
    !-------------------------------------------------------------------
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
           DensityA = C0*DensityA_T + (1.0d0-C0)*DensityA
           DensityB = C0*DensityB_T + (1.0d0-C0)*DensityB
       end if
       
       !------------------------ G matrix ------------------------------
       call UHF_Const_G(N,Nmem,Neri,Two_Electron,EI,DensityA,DensityB,Coul,ExcA,ExcB,t2,d2,CE)
       !----------------------------------------------------------------
        
       !------------------- rho,Func_E,Vrho,Vsigma ---------------------
       call rho_Func_Open(Node1,M,DensityA,DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,method, &
                          gridx(1:Node1,1:M),gridy(1:Node1,1:M),gridz(1:Node1,1:M), &
                          end_pot,Rho_gridA(1:Node1,1:M,1:3),Rho_gridB(1:Node1,1:M,1:3), &
                          Func_E(1:Node1,1:M),VrhoA(1:Node1,1:M),VrhoB(1:Node1,1:M), &
                          VsigmaA(1:Node1,1:M),VsigmaB(1:Node1,1:M),VsigmaAB(1:Node1,1:M)) 

       !------------------------- XC Matrix ----------------------------       
       XCA = 0.0d0 ; XCB = 0.0d0
       if( method == "XAlpha" .or. method == "xalpha" &
           .or. method == "LSDA" .or. method == "lsda" ) then
      
         do i = 1,M
            do j = 1,end_pot(i)    
              x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)  
              !---------------------------------------------------------
              do i1 = 1,N              
                 Phi(i1) = 0.0d0
                 D_Phi(:,i1) = 0.0d0
                 do k1 = 1,bas_num(i1)
                      a = a_Basis(k1,i1,NAtom(i1)) ; Ra = Rt(:,i1) 
                      La = l_Basis(:,i1)
                      call Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,"N")
                      Phi(i1) = Phi(i1) + ATemp0*c_Basis(k1,i1,NAtom(i1))
                 end do
               end do             
               !--------------------------------------------------------
               do k = 1,N
                  do l = k,N
                  
                      XCA(l,k) =  XCA(l,k) + Phi(k)*Phi(l)*VrhoA(j,i)*W_i(j,i)
                      XCB(l,k) =  XCB(l,k) + Phi(k)*Phi(l)*VrhoB(j,i)*W_i(j,i)

                  end do       
               end do 
               !--------------------------------------------------------
            end do
         end do 
       
       else
      
         do i = 1,M
            do j = 1,end_pot(i)  
              x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)  
              !---------------------------------------------------------
              do i1 = 1,N              
                 Phi(i1) = 0.0d0
                 D_Phi(:,i1) = 0.0d0
                 do k1 = 1,bas_num(i1)
                      a = a_Basis(k1,i1,NAtom(i1)) ; Ra = Rt(:,i1) 
                      La = l_Basis(:,i1)
                      call Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,"Y")
                      Phi(i1) = Phi(i1) + ATemp0*c_Basis(k1,i1,NAtom(i1))
                      D_Phi(:,i1) = D_Phi(:,i1) + A3Temp0*c_Basis(k1,i1,NAtom(i1))
                 end do
               end do             
               !---------------------------------------------------------    
               do k = 1,N
                  do l = k,N
                      
                      Temp1 = Phi(k)*Phi(l)*VrhoA(j,i) 
                      Temp2 = sum( ( 2.0d0*VsigmaA(j,i)*Rho_gridA(j,i,:) + VsigmaAB(j,i)*Rho_gridB(j,i,:) )*&
                              ( Phi(k)*D_Phi(:,l) + Phi(l)*D_Phi(:,k) ) )
                      XCA(l,k) =  XCA(l,k) + ( Temp1 + Temp2 )*W_i(j,i)
                      
                      Temp1 = Phi(k)*Phi(l)*VrhoB(j,i) 
                      Temp2 = sum( ( 2.0d0*VsigmaB(j,i)*Rho_gridB(j,i,:) + VsigmaAB(j,i)*Rho_gridA(j,i,:) )*&
                              ( Phi(k)*D_Phi(:,l) + Phi(l)*D_Phi(:,k) ) )
                      XCB(l,k) =  XCB(l,k) + ( Temp1 + Temp2 )*W_i(j,i)
                                  
                  end do       
               end do
               !-------------------------------------------------------- 
            end do
         end do  
       
       end if  
       
       do k = 1,N
          do l = k,N
              XCA(k,l) = XCA(l,k) 
              XCB(k,l) = XCB(l,k)  
          end do
       end do                 
       !----------------------------------------------------------------
       E_Exc = 0.0d0 ; w0 = 0.0d0
       if( method == "rob3lyp" .or. method == "ropbe1pbe" .or. method == "rob3pw91" ) then
           
          if( method == "rob3lyp" ) w0 = 0.2d0
          if( method == "rob3pw91" ) w0 = 0.2d0
          if( method == "ropbe1pbe" ) w0 = 0.25d0
         
          XCA = XCA + w0*ExcA
          XCB = XCB + w0*ExcB 
          
          Temp1 = 0.0d0
          do i = 1,N
              Temp2 = 0.0d0
              do j = 1,N
                   Temp2 = Temp2 + DensityA(j,i)*ExcA(j,i)+DensityB(j,i)*ExcB(j,i)
              end do
              Temp1 = Temp1 + Temp2 
          end do     
          E_Exc = Temp1*0.5d0
               
       end if
       !----------------------------------------------------------------
       
       idi = mod(ii,7)
       if( mod(ii,7) == 0 ) idi = 7
                                           
       FockA = H_core + Coul + XCA
       FockB = H_core + Coul + XCB        !Fock Matrix
       
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
           if( ii <= 7 ) then
                call UC_Coff(idi,N,ErA(:,:,1:idi),ErB(:,:,1:idi),CCoff(1:idi+1))
                Fock(:,:,0) = 0.0d0
                do i = 1,idi
                   Fock(:,:,0) = Fock(:,:,0) + Fock(:,:,i)*CCoff(i)
                end do 
           else
                call UC_Coff(7,N,ErA(:,:,1:7),ErB(:,:,1:7),CCoff(1:8))
                Fock(:,:,0) = 0.0d0
                do i = 1,7
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
                Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*(H_core(j,i) + 0.5d0*Coul(j,i))
           end do
           Temp1 = Temp1 + Temp2 
       end do
       
       call Int_Func_E_XC(Node1,M,Func_E(1:Node1,1:M),W_i(1:Node1,1:M),end_pot,E_XC)
             
       Ee = Temp1 + E_XC + w0*E_Exc           !Calculate Energy E
       Delta_E = dabs( Ee - ETemp )
       ETemp = Ee
       
       write(*,"(A14,ES16.8)") "E(DFT)= ",Ee + E_Nuc 
       write(*,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(*,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(100,"(A14,ES16.8)") "E(DFT)= ",Ee + E_Nuc 
       write(100,"(A14,ES14.6)") "DIIS_ER=",max(emaxA,emaxB)
       write(100,"(A13,ES15.6,A10,ES16.6)") "DDen  =",max(Delta_DenA,Delta_DenB),"DE =",Delta_E
       write(*,*) 
       write(100,*)
       if( Delta_E < 10.0d0**(-Conv) .and. Delta_DenA < 1.0d-4 .and. Delta_DenB < 1.0d-4) exit
       !Judge whether the SCF procedure is converged.
       
       if( Delta_E < 10.0d0**(-5.0d0) .and. jud_node) then 
           jud_node = .false.
           Node1 = Node
           write(*,*) "     Initial convergence to 1.0D-05 achieved. "
           write(*,*) "     Increase integral accuracy."
           write(*,*)
           write(100,*) "     Initial convergence to 1.0D-05 achieved. "
           write(100,*) "     Increase integral accuracy."
           write(100,*)
           call Int_Dep(Node1,M,Rn,MAtom,gridx(1:Node1,1:M),gridy(1:Node1,1:M), &
                gridz(1:Node1,1:M), W_i(1:Node1,1:M),end_pot)
       end if
              
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
    
end subroutine DFT_RO
!=======================================================================




