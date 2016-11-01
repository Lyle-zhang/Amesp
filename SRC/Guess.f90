subroutine AtDen(M,N,Nmem,Neri,MAtom,NAtom,Overlap,Kinetic,H_core,Rn,Rt,a_Basis,&
                 c_Basis,l_Basis,bas_num,Two_Electron,EI,bas_set,t2,d2,MoCu)
!-----------------------------------------------------------------------
    ! Use natural orbitals of a diagonal density matrix constructed 
    ! using atomic orbitals and atomic occupation number.
!-----------------------------------------------------------------------
    implicit none
    integer ii,jj,i1,i2,j1,j2,j3,i,j,k,l,M,N,Nmem,Neri,EI(N,N)
    integer ,allocatable :: EI0(:,:)
    integer N_bas(36),N_occA(36),N_occB(36),MAtom(M),l_Basis(3,N)
    integer bas_num(N),NAtom(N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36),Rn(3,M),Rt(3,N)
    real(kind=8) MoCu(N,N),Overlap_Temp(N,N),Sqrt_S(N,N),Fock(N,N)
    real(kind=8) Overlap(N,N),H_Core(N,N),Kinetic(N,N),HCore(42,42)
    real(kind=8) N_E_Attract(42,42),Density(N,N),En(N),G(N,N),EnA(42),EnB(42)
    real(kind=8) DensityA(N,N),DensityB(N,N),GA(N,N),GB(N,N)
    real(kind=8) FockA(42,42),FockB(42,42),MoCuA(42,42),MoCuB(42,42)
    real(kind=8) ,allocatable :: Two_Ele(:,:)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) Temp1,Temp2,Delta_E,Ee,ETemp,EeA,EeB
    character(len=20) bas_set
    character(len=10) t2
    character(len=8) d2

    if( bas_set == "sto-3g" ) then    
         N_bas = (/1,1,5,5,5,5,5,5,5,5,9,9,9,9,9,9,9,9,13,13,19,19,19,19,19, &
                 19,19,19,19,19,19,19,19,19,19,19/)
    else if( bas_set == "3-21g"  ) then   
         N_bas = (/2,2,9,9,9,9,9,9,9,9,13,13,13,13,13,13,13,13,17,17,29,29, &
                 29,29,29,29,29,29,29,29,23,23,23,23,23,23/)
    else if( bas_set == "6-31g" ) then
         N_bas = (/2,2,9,9,9,9,9,9,9,9,13,13,13,13,13,13,13,13, &
                  17,17,29,29,29,29,29,29,29,29,29,29,24,24,24,24,24,24/)
    else if( bas_set == "6-31g*" ) then
         N_bas(1:20) = (/2,2,15,15,15,15,15,15,15,15,19,19,19, &
                       9,19,19,19,19,23,23/)
    else if( bas_set == "6-31g**" ) then
         N_bas(1:20) = (/5,5,15,15,15,15,15,15,15,15,19,19,19, &
                        19,19,19,19,19,23,23/)
    else if( bas_set == "6-311g"  ) then   
         N_bas(1:20) = (/3,3,13,13,13,13,13,13,13,13,21,21,21,21,21,21,21,21,35,35/)
    end if
    N_occA(1:36) = (/1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12, &
                     13,13,14,14,15,15,16,16,17,17,18,18/)
    N_occB(1:36) = (/0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12, &
                     12,13,13,14,14,15,15,16,16,17,17,18/)  
    
    Density = 0.0d0 ; DensityA = 0.0d0 ; DensityB = 0.0d0
    
    i1 = 0
    do ii = 1,M
    
          j1 = N_bas(MAtom(ii))
          j2 = N_occA(MAtom(ii))
          j3 = N_occB(MAtom(ii))
          i2 = i1 + j1
          
          allocate(Two_Ele(j1*(j1+1)/2,j1*(j1+1)/2),EI0(j1,j1))
          
          !============================RHF==============================
          if( mod(MAtom(ii),2) == 0 ) then
              Overlap_Temp(1:j1,1:j1) = Overlap(i1+1:i2,i1+1:i2)
              call Mat_Sqrt(Overlap_Temp(1:j1,1:j1),Sqrt_S(1:j1,1:j1),j1)
              
              call N_E_Attractive_Int(N_E_Attract(1:j1,1:j1),1,j1,a_Basis(:,i1+1:i2,:), &
                                      c_Basis(:,i1+1:i2,:),l_Basis(:,i1+1:i2),bas_num(i1+1:i2), &
                                      Rn(:,ii),Rt(:,i1+1:i2),MAtom(ii),NAtom(i1+1:i2))
              call Two_Electron_Int(Two_Ele,EI0,j1,j1,j1,a_Basis(:,i1+1:i2,:),c_Basis(:,i1+1:i2,:), &
                                      l_Basis(:,i1+1:i2),bas_num(i1+1:i2),Rt(:,i1+1:i2),NAtom(i1+1:i2),t2,d2)
                                      
              HCore(1:j1,1:j1) = Kinetic(i1+1:i2,i1+1:i2) + N_E_Attract(1:j1,1:j1)
              
              call Main_Calculate_S(Sqrt_S(1:j1,1:j1),HCore(1:j1,1:j1),MoCu(1:j1,1:j1),En(1:j1),j1)              
              
              ETemp = 1.0d7
              do jj = 1,128
              !---------------------------------------------------------
                   Density(i1+1:i2,i1+1:i2) = matmul(MoCu(1:j1,1:j2),transpose(MoCu(1:j1,1:j2)))
                   
                   do i = 1,j1
                      do j = i,j1
                         Temp1 = 0.0d0
                         do k = 1,j1
                           Temp2 = 0.0d0
                           do l = 1,j1
                             Temp2 = Temp2 + (2.0d0*Two_Ele(EI0(l,k),EI0(j,i))- &
                                     Two_Ele(EI0(j,k),EI0(l,i)))*Density(k+i1,l+i1)
                           end do
                           Temp1 = Temp1 + Temp2
                         end do
                         G(j,i) = Temp1      !Calculate G(j,i)
                         G(i,j) = Temp1
                      end do
                  end do
              
                  Fock(1:j1,1:j1) = HCore(1:j1,1:j1) + G(1:j1,1:j1)
                  call Main_Calculate_S(Sqrt_S(1:j1,1:j1),Fock(1:j1,1:j1), &
                       MoCu(1:j1,1:j1),En(1:j1),j1) !Solve FC=SCE
                  
                  Temp1 = 0.0d0
                  do i = i1+1,i2
                       Temp2 = 0.0d0
                       do j = i1+1,i2
                            Temp2 = Temp2 + Density(j,i)*(Hcore(j-i1,i-i1) + Fock(j-i1,i-i1))
                       end do
                       Temp1 = Temp1 + Temp2 
                  end do
                  Ee = Temp1                 !Calculate Energy E
                  Delta_E = dabs( Ee - ETemp )
                  if( Delta_E < 1.0d-8 ) exit 
                  !Judge whether the SCF procedure is converged.
                  ETemp = Ee
              !---------------------------------------------------------
              end do
              Density(i1+1:i2,i1+1:i2) = matmul(MoCu(1:j1,1:j2),transpose(MoCu(1:j1,1:j2)))          
          end if
          !============================UHF==============================
          if( mod(MAtom(ii),2) == 1 ) then
              Overlap_Temp(1:j1,1:j1) = Overlap(i1+1:i2,i1+1:i2)
              call Mat_Sqrt(Overlap_Temp(1:j1,1:j1),Sqrt_S(1:j1,1:j1),j1)
                            
              call N_E_Attractive_Int(N_E_Attract(1:j1,1:j1),1,j1,a_Basis(:,i1+1:i2,:), &
                                      c_Basis(:,i1+1:i2,:),l_Basis(:,i1+1:i2),bas_num(i1+1:i2), &
                                      Rn(:,ii),Rt(:,i1+1:i2),MAtom(ii),NAtom(i1+1:i2))
              call Two_Electron_Int(Two_Ele,EI0,j1,j1,j1,a_Basis(:,i1+1:i2,:),c_Basis(:,i1+1:i2,:), &
                                      l_Basis(:,i1+1:i2),bas_num(i1+1:i2),Rt(:,i1+1:i2),NAtom(i1+1:i2),t2,d2)
                                      
              HCore(1:j1,1:j1) = Kinetic(i1+1:i2,i1+1:i2) + N_E_Attract(1:j1,1:j1)
              
              call Main_Calculate_S(Sqrt_S(1:j1,1:j1),HCore(1:j1,1:j1),MoCuA(1:j1,1:j1),EnA(1:j1),j1)   
              MoCuB(1:j1,1:j1) = MoCuA(1:j1,1:j1)              
              
              ETemp = 1.0d7
              do jj = 1,128
              
                  DensityA(i1+1:i2,i1+1:i2) = matmul(MoCuA(1:j1,1:j2),transpose(MoCuA(1:j1,1:j2)))  
                  DensityB(i1+1:i2,i1+1:i2) = matmul(MoCuB(1:j1,1:j3),transpose(MoCuB(1:j1,1:j3))) 
                   
                   do i = 1,j1
                      do j = i,j1
                         Temp1 = 0.0d0
                         do k = 1,j1
                           Temp2 = 0.0d0
                           do l = 1,j1
                                Temp2 = Temp2 + (DensityA(k+i1,l+i1)+DensityB(k+i1,l+i1))* &
                                        Two_Ele(EI0(l,k),EI0(j,i)) - &
                                        DensityA(k+i1,l+i1)*Two_Ele(EI0(j,k),EI0(l,i))
                           end do
                           Temp1 = Temp1 + Temp2
                         end do
                         GA(j,i) = Temp1      !Calculate GA(j,i)
                         GA(i,j) = Temp1
                      end do
                  end do
                   
                   do i = 1,j1
                      do j = i,j1
                         Temp1 = 0.0d0
                         do k = 1,j1
                           Temp2 = 0.0d0
                           do l = 1,j1
                                Temp2 = Temp2 + (DensityA(k+i1,l+i1)+DensityB(k+i1,l+i1))* &
                                        Two_Ele(EI0(l,k),EI0(j,i)) - &
                                        DensityB(k+i1,l+i1)*Two_Ele(EI0(j,k),EI0(l,i))
                           end do
                           Temp1 = Temp1 + Temp2
                         end do
                         GB(j,i) = Temp1      !Calculate GB(j,i)
                         GB(i,j) = Temp1
                      end do
                  end do
                   
                  FockA(1:j1,1:j1) = HCore(1:j1,1:j1) + GA(1:j1,1:j1)
                  FockB(1:j1,1:j1) = HCore(1:j1,1:j1) + GB(1:j1,1:j1)
  
                  call Main_Calculate_S(Sqrt_S(1:j1,1:j1),FockA(1:j1,1:j1), &
                       MoCuA(1:j1,1:j1),EnA(1:j1),j1) !Solve FC=SCE
                  call Main_Calculate_S(Sqrt_S(1:j1,1:j1),FockB(1:j1,1:j1), &
                       MoCuB(1:j1,1:j1),EnB(1:j1),j1) !Solve FC=SCE
                   
                  Temp1 = 0.0d0
                  do i = i1+1,i2
                      Temp2 = 0.0d0
                      do j = i1+1,i2
                          Temp2 = Temp2 + (DensityA(j,i)+DensityB(j,i))*Hcore(j-i1,i-i1)
                      end do
                      Temp1 = Temp1 + Temp2 
                  end do
       
                  EeA = 0.0d0
                  do i = 1,j2
                     EeA = EeA + EnA(i)
                  end do
       
                  EeB = 0.0d0
                  do i = 1,j3
                     EeB = EeB + EnB(i)
                  end do
       
                  Ee = 0.5d0*( Temp1 + EeA + EeB )               !Calculate Energy E
                  Delta_E = dabs( Ee - ETemp )           
                  ETemp = Ee
                  if( Delta_E < 1.0d-8 ) exit 
                  !Judge whether the SCF procedure is converged. 
                                                               
              end do          
              DensityA(i1+1:i2,i1+1:i2) = matmul(MoCuA(1:j1,1:j2),transpose(MoCuA(1:j1,1:j2)))  
              DensityB(i1+1:i2,i1+1:i2) = matmul(MoCuB(1:j1,1:j3),transpose(MoCuB(1:j1,1:j3)))
              Density(i1+1:i2,i1+1:i2) = 0.5d0*(DensityA(i1+1:i2,i1+1:i2)+DensityB(i1+1:i2,i1+1:i2))              
          end if
          !=============================================================
          
          i1 = i2 
          deallocate(Two_Ele,EI0)    
          
    end do
!-----------------------------------------------------------------------

    call RHF_Const_G(N,Nmem,Neri,Two_Electron,EI,Density,G,t2,d2)
    
    Fock = H_core + G
    
    Overlap_Temp = Overlap
    call Mat_Sqrt(Overlap_Temp,Sqrt_S,N)
    
    call Main_Calculate_S(Sqrt_S,Fock,MoCu,En,N) !Solve FC=SCE
     
end subroutine AtDen
!=======================================================================
