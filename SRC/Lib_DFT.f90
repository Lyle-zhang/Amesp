subroutine Int_Dep(Node,M,Rn,MAtom,gridx,gridy,gridz,W_i,end_pot)
    ! Calculate the grids and weight for DFT.
    implicit none
    integer Node,Nsph1,Nsph2,Nrad1,Nrad2,Nrad3
    integer M,i,j,k,ii,jj,NN
    integer rad_pot(M),sph_pot(M),r_cut(M),end_pot(M)
    integer Bond(M,M),MAtom(M)
    real(kind=8),parameter::Pi = 3.14159265358979d0
    real(kind=8) chi,uij,aij,nu_ij
    real(kind=8) rcov(20),threshold
    real(kind=8) R,x,y,z,rj,xj,wj
    real(kind=8) Rn(3,M),Distance(M,M),r_i,r_j,miu_ij
    real(kind=8) S_u(M,M),Temp,W_Temp(M),W_i0
    real(kind=8),allocatable:: X1(:,:),X2(:,:)
    real(kind=8),allocatable:: pot_x(:),pot_y(:),pot_z(:),pot_w(:)
    real(kind=8) grid0x(Node,M),grid0y(Node,M),grid0z(Node,M),gridw(Node,M)
    real(kind=8) W_i(Node,M),gridx(Node,M),gridy(Node,M),gridz(Node,M)

    if( Node == 2960 ) then 
        Nsph1 = 50
        Nsph2 = 74
        Nrad1 = 25   
        Nrad2 = 35  
        Nrad3 = 40 
        threshold = 1.0d-10  
    end if
       
    if( Node == 8760 ) then 
        Nsph1 = 110
        Nsph2 = 146
        Nrad1 = 30   
        Nrad2 = 45  
        Nrad3 = 60 
        threshold = 1.0d-11 
    end if
    
    if( Node == 22650 ) then 
        Nsph1 = 266
        Nsph2 = 302
        Nrad1 = 35   
        Nrad2 = 60  
        Nrad3 = 75 
        threshold = 1.0d-12
    end if
    
    allocate(X1(4,Nsph2),X2(4,Nsph2),pot_x(Nsph2),pot_y(Nsph2),pot_z(Nsph2),pot_w(Nsph2))
    
    !rcov = (/0.32d0,0.46d0 &
            !,1.33d0,1.02d0,0.85d0,0.75d0,0.71d0,0.63d0,0.64d0,0.96d0 &
            !,1.55d0,1.39d0,1.26d0,1.16d0,1.11d0,1.03d0,0.99d0,0.96d0 &
            !,1.96d0,1.71d0/)/0.5291772083d0         !Covalent Radius

    rcov = (/0.5291772083d0,0.31d0 &
            ,1.67d0,1.12d0,0.87d0,0.67d0,0.56d0,0.48d0,0.42d0,0.38d0 &
            ,1.90d0,1.45d0,1.18d0,1.11d0,0.98d0,0.88d0,0.79d0,0.71d0 &
            ,2.43d0,1.94d0/)/0.5291772083d0         !Covalent Radius
               
    Bond = 0
    Distance = 0.0d0
    do i = 1,M
       do j = i+1,M
          Distance(i,j) = dsqrt((Rn(1,i)-Rn(1,j))**2+(Rn(2,i)-Rn(2,j))**2+(Rn(3,i)-Rn(3,j))**2)
          Distance(j,i) = Distance(i,j)
          if( distance(i,j) < 1.15d0*(rcov(MAtom(i))+rcov(MAtom(j))) ) then
              Bond(j,i) = 1
              Bond(i,j) = 1 
             !distance < 1.15*(rcov(i)+rcov(j)) , rcov(H)=0.35/0.5291772083
          end if
       end do
    end do
    
    do i = 1,M
    
        if ( MAtom(i) <= 2 ) rad_pot(i) =  Nrad1
        if ( MAtom(i) <= 10 .and. MAtom(i) > 2 ) rad_pot(i) = Nrad2
        if ( MAtom(i) > 10 ) rad_pot(i) = Nrad3 !75 for Gaussian.
        
        if (sum(Bond(:,i)) == 1) Then !302 for Gaussian.
           sph_pot(i) = Nsph1
           
           if ( Nsph1 == 50 ) then
                call LD0050(X1(1,1:Nsph1),X1(2,1:Nsph1),X1(3,1:Nsph1),X1(4,1:Nsph1),NN)
           else if( Nsph1 == 110 ) then
                call LD0110(X1(1,1:Nsph1),X1(2,1:Nsph1),X1(3,1:Nsph1),X1(4,1:Nsph1),NN)
           else if( Nsph1 == 266 ) then
                call LD0266(X1(1,1:Nsph1),X1(2,1:Nsph1),X1(3,1:Nsph1),X1(4,1:Nsph1),NN)
           end if
           
           pot_x(1:sph_pot(i)) = X1(1,1:Nsph1) ; pot_y(1:sph_pot(i)) = X1(2,1:Nsph1)
           pot_z(1:sph_pot(i)) = X1(3,1:Nsph1) ; pot_w(1:sph_pot(i)) = X1(4,1:Nsph1)
        else
           sph_pot(i) = Nsph2

           if ( Nsph2 == 74 ) then
                call LD0074(X2(1,1:Nsph2),X2(2,1:Nsph2),X2(3,1:Nsph2),X2(4,1:Nsph2),NN)           
           else if( Nsph2 == 146 ) then
                call LD0146(X2(1,1:Nsph2),X2(2,1:Nsph2),X2(3,1:Nsph2),X2(4,1:Nsph2),NN)
           else if( Nsph2 == 302 ) then
                call LD0302(X2(1,1:Nsph2),X2(2,1:Nsph2),X2(3,1:Nsph2),X2(4,1:Nsph2),NN)
           end if
           
           pot_x(1:sph_pot(i)) = X2(1,1:Nsph2) ; pot_y(1:sph_pot(i)) = X2(2,1:Nsph2)
           pot_z(1:sph_pot(i)) = X2(3,1:Nsph2) ; pot_w(1:sph_pot(i)) = X2(4,1:Nsph2) 
        end if   
   
        R = rcov(MAtom(i))/2.0d0              !(1+x_i)/(1-x_i)*R
        if( MAtom(i) == 1 ) R = rcov(MAtom(i))
        r_cut(i) = 0
        
        !================= gridx,gridy,gridz,gridw =====================
        do j = 1,rad_pot(i)
           !------------------------------------------------------ 
           ! Gauss-Chebyshev
           xj = dcos(j*Pi/(rad_pot(i)+1.0d0))
           rj = R*(1.0d0+xj)/(1.0d0-xj)
           wj = 2.0d0*Pi/(rad_pot(i)+1.0d0)*R**3*(1.0d0+xj)**2.5d0/(1.0d0-xj)**3.5d0
           !---------------------------------------------------------------
           ! Euler-Maclaurin
           !xj = j/(rad_pot(i)+1.0d0)
           !rj = R*xj**2/(1.0d0-xj)**2
           !wj = 2.0d0*R**3/(rad_pot(i)+1.0d0)*xj**5.0d0/(1.0d0-xj)**7
           !------------------------------------------------------ 
           grid0x((j-1)*sph_pot(i)+1:j*sph_pot(i),i) = rj*pot_x
           grid0y((j-1)*sph_pot(i)+1:j*sph_pot(i),i) = rj*pot_y 
           grid0z((j-1)*sph_pot(i)+1:j*sph_pot(i),i) = rj*pot_z 
           gridw((j-1)*sph_pot(i)+1:j*sph_pot(i),i) = wj*pot_w*4.0d0*Pi
           if( r_cut(i) == 0 .and. rj < 18.897261339d0 ) r_cut(i) = j
        end do
        
        grid0x(:,i) = grid0x(:,i) + Rn(1,i)
        grid0y(:,i) = grid0y(:,i) + Rn(2,i)
        grid0z(:,i) = grid0z(:,i) + Rn(3,i)
        !===============================================================
        k = 0
        do j = 1+r_cut(i)*sph_pot(i),rad_pot(i)*sph_pot(i)
           
           x = grid0x(j,i) ; y = grid0y(j,i) ; z = grid0z(j,i)       
           !========================== W_i =============================
           S_u = 1.0d0
           do ii = 1,M
              r_i = dsqrt((x-Rn(1,ii))**2+(y-Rn(2,ii))**2+(z-Rn(3,ii))**2)
              do jj = 1,M              
                 if(ii == jj) cycle
                 r_j = dsqrt((x-Rn(1,jj))**2+(y-Rn(2,jj))**2+(z-Rn(3,jj))**2)
                 miu_ij = ( r_i - r_j )/Distance(jj,ii) !μ(i,j)
                 
                 !Change μ(i,j) to ν(i,j)
                 chi = rcov(MAtom(ii))/rcov(MAtom(jj))
                 if( chi < 1.0d0+1d-6 .and. chi > 1.0d0-1d-6) then
                    nu_ij = miu_ij
                 else
                    uij = (chi-1.0d0)/(chi+1.0d0)
                    aij = uij/(uij**2-1.0d0)
                    if (aij > 0.5D0) aij = 0.5D0
                    if (aij < -0.5D0) aij = -0.5D0
                    nu_ij = miu_ij + aij*(1.0d0 - miu_ij**2) !ν(i,j)
                 end if
                                  
                 Temp = 1.5d0*nu_ij - 0.5d0*nu_ij**3
                 Temp = 1.5d0*Temp - 0.5d0*Temp**3
                 Temp = 1.5d0*Temp - 0.5d0*Temp**3
                 S_u(jj,ii) = 0.5d0*( 1.0d0 - Temp ) !s(μ(i,j))=0.5(1-p(p(p(ν(i,j)))))                
              end do
           end do
           
           W_Temp = 1.0d0
           do ii = 1,M
              W_Temp = W_Temp*S_u(ii,:)
           end do
           W_i0 = gridw(j,i)*W_Temp(i)/sum(W_Temp)

           if ( dabs(W_i0) > threshold ) then
              k = k + 1
              W_i(k,i) = W_i0
              gridx(k,i) = x
              gridy(k,i) = y
              gridz(k,i) = z
           end if
           
        end do   !do j = 1+r_cut*sph_pot,rad_pot*sph_pot
        end_pot(i) = k
    end do   !do i = 1,M

end subroutine Int_Dep
!=======================================================================
subroutine rho_Func_Close(Node,M,Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,method, &
               gridx,gridy,gridz,end_pot,Rho_grid,Func_E,Vrho,Vsigma)
    ! rho,Func_E,Vrho,Vsigma                  
    implicit none 
    integer Node
    integer M,N,i,j
    integer end_pot(M),NAtom(N),bas_num(N),l_Basis(3,N)
    real(kind=8) x,y,z
    real(kind=8) c_Basis(3,N,2),a_Basis(3,N,2),Rt(3,N),Density(N,N)
    real(kind=8) gridx(Node,M),gridy(Node,M),gridz(Node,M)
    real(kind=8) Func_E(Node,M),Vrho(Node,M),Vsigma(Node,M),rho,Rho_grid(Node,M,3),sigma
    character(len=20) method
    
    !--------------------- rho,Func_E,Vrho,Vsigma ----------------------
    if( method == "XAlpha" .or. method == "xalpha" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"N")
             call XAlpha_Close(rho,Func_E(j,i),Vrho(j,i))  
              
          end do
       end do
       Vsigma = 0.0d0
    end if
    !----------------------------------------------------------------  
    if( method == "LSDA" .or. method == "lsda" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"N")
             call LSDA_Close(rho,Func_E(j,i),Vrho(j,i))  
              
          end do
       end do
       Vsigma = 0.0d0
    end if
    !----------------------------------------------------------------  
    if( method == "PBEPBE" .or. method == "pbepbe" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigma = Rho_grid(j,i,1)**2 + Rho_grid(j,i,2)**2 + Rho_grid(j,i,3)**2
             call pbe_Close(rho,sigma,Func_E(j,i),Vrho(j,i),Vsigma(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------  
    if( method == "PW91PW91" .or. method == "pw91pw91" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigma = Rho_grid(j,i,1)**2 + Rho_grid(j,i,2)**2 + Rho_grid(j,i,3)**2
             call pw91_Close(rho,sigma,Func_E(j,i),Vrho(j,i),Vsigma(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "b3lyp" .or. method == "B3LYP" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigma = Rho_grid(j,i,1)**2 + Rho_grid(j,i,2)**2 + Rho_grid(j,i,3)**2
             call B3LYP_Close(rho,sigma,Func_E(j,i),Vrho(j,i),Vsigma(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "blyp" .or. method == "BLYP" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigma = Rho_grid(j,i,1)**2 + Rho_grid(j,i,2)**2 + Rho_grid(j,i,3)**2
             call blyp_Close(rho,sigma,Func_E(j,i),Vrho(j,i),Vsigma(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "pbe1pbe" .or. method == "PBE1PBE" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigma = Rho_grid(j,i,1)**2 + Rho_grid(j,i,2)**2 + Rho_grid(j,i,3)**2
             call pbe0_Close(rho,sigma,Func_E(j,i),Vrho(j,i),Vsigma(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "b3pw91" .or. method == "B3PW91" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rho,Rho_grid(j,i,:),Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigma = Rho_grid(j,i,1)**2 + Rho_grid(j,i,2)**2 + Rho_grid(j,i,3)**2
             call b3pw91_Close(rho,sigma,Func_E(j,i),Vrho(j,i),Vsigma(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
         
end subroutine rho_Func_Close
!=======================================================================
subroutine rho_Func_Open(Node,M,DensityA,DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,method, &
                         gridx,gridy,gridz,end_pot,Rho_gridA,Rho_gridB, &
                         Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)
    ! rho,Func_E,Vrho,Vsigma                      
    implicit none 
    integer Node
    integer M,N,i,j
    integer end_pot(M),NAtom(N),bas_num(N),l_Basis(3,N)
    real(kind=8) x,y,z,Rho_gridA(Node,M,3),Rho_gridB(Node,M,3)
    real(kind=8) VsigmaAB(Node,M),VsigmaA(Node,M),VsigmaB(Node,M)
    real(kind=8) c_Basis(3,N,2),a_Basis(3,N,2),Rt(3,N),DensityA(N,N),DensityB(N,N)
    real(kind=8) gridx(Node,M),gridy(Node,M),gridz(Node,M),sigmaAB
    real(kind=8) Func_E(Node,M),VrhoA(Node,M),rhoA,sigmaA,VrhoB(Node,M),rhoB,sigmaB
    character(len=20) method
    
    !--------------------- rho,Func_E,Vrho,Vsigma ----------------------
    if( method == "XAlpha" .or. method == "xalpha" .or. method == "roxalpha") then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"N")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"N")
             call XAlpha_Open(rhoA,rhoB,Func_E(j,i),VrhoA(j,i),VrhoB(j,i)) 
              
          end do
       end do
       VsigmaA = 0.0d0 ; VsigmaB = 0.0d0
    end if
    !---------------------------------------------------------------- 
    if( method == "LSDA" .or. method == "lsda" .or. method == "rolsda") then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"N")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"N")
             call LSDA_Open(rhoA,rhoB,Func_E(j,i),VrhoA(j,i),VrhoB(j,i)) 
              
          end do
       end do
       VsigmaA = 0.0d0 ; VsigmaB = 0.0d0
    end if
    !---------------------------------------------------------------- 
    if( method == "PBEPBE" .or. method == "pbepbe" .or. method == "ropbepbe") then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigmaA = Rho_gridA(j,i,1)**2 + Rho_gridA(j,i,2)**2 + Rho_gridA(j,i,3)**2
             sigmaB = Rho_gridB(j,i,1)**2 + Rho_gridB(j,i,2)**2 + Rho_gridB(j,i,3)**2
             sigmaAB = Rho_gridB(j,i,1)*Rho_gridA(j,i,1) + Rho_gridB(j,i,2)*Rho_gridA(j,i,2) + &
                       Rho_gridB(j,i,3)*Rho_gridA(j,i,3)
             call pbe_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E(j,i),VrhoA(j,i), &
                           VrhoB(j,i),VsigmaA(j,i),VsigmaB(j,i),VsigmaAB(j,i))
              
          end do
       end do
    end if
    !---------------------------------------------------------------- 
    if( method == "PW91PW91" .or. method == "pw91pw91" .or. method == "ropw91pw91") then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigmaA = Rho_gridA(j,i,1)**2 + Rho_gridA(j,i,2)**2 + Rho_gridA(j,i,3)**2
             sigmaB = Rho_gridB(j,i,1)**2 + Rho_gridB(j,i,2)**2 + Rho_gridB(j,i,3)**2
             sigmaAB = Rho_gridB(j,i,1)*Rho_gridA(j,i,1) + Rho_gridB(j,i,2)*Rho_gridA(j,i,2) + &
                       Rho_gridB(j,i,3)*Rho_gridA(j,i,3)
             call pw91_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E(j,i),VrhoA(j,i), &
                           VrhoB(j,i),VsigmaA(j,i),VsigmaB(j,i),VsigmaAB(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "b3lyp" .or. method == "B3LYP" .or. method == "rob3lyp" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigmaA = Rho_gridA(j,i,1)**2 + Rho_gridA(j,i,2)**2 + Rho_gridA(j,i,3)**2
             sigmaB = Rho_gridB(j,i,1)**2 + Rho_gridB(j,i,2)**2 + Rho_gridB(j,i,3)**2
             sigmaAB = Rho_gridB(j,i,1)*Rho_gridA(j,i,1) + Rho_gridB(j,i,2)*Rho_gridA(j,i,2) + &
                       Rho_gridB(j,i,3)*Rho_gridA(j,i,3)
             call B3LYP_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E(j,i),VrhoA(j,i), &
                           VrhoB(j,i),VsigmaA(j,i),VsigmaB(j,i),VsigmaAB(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "blyp" .or. method == "BLYP" .or. method == "roblyp") then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigmaA = Rho_gridA(j,i,1)**2 + Rho_gridA(j,i,2)**2 + Rho_gridA(j,i,3)**2
             sigmaB = Rho_gridB(j,i,1)**2 + Rho_gridB(j,i,2)**2 + Rho_gridB(j,i,3)**2
             sigmaAB = Rho_gridB(j,i,1)*Rho_gridA(j,i,1) + Rho_gridB(j,i,2)*Rho_gridA(j,i,2) + &
                       Rho_gridB(j,i,3)*Rho_gridA(j,i,3)
             call blyp_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E(j,i),VrhoA(j,i), &
                           VrhoB(j,i),VsigmaA(j,i),VsigmaB(j,i),VsigmaAB(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "pbe1pbe" .or. method == "PBE1PBE" .or. method == "ropbe1pbe") then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigmaA = Rho_gridA(j,i,1)**2 + Rho_gridA(j,i,2)**2 + Rho_gridA(j,i,3)**2
             sigmaB = Rho_gridB(j,i,1)**2 + Rho_gridB(j,i,2)**2 + Rho_gridB(j,i,3)**2
             sigmaAB = Rho_gridB(j,i,1)*Rho_gridA(j,i,1) + Rho_gridB(j,i,2)*Rho_gridA(j,i,2) + &
                       Rho_gridB(j,i,3)*Rho_gridA(j,i,3)
             call pbe0_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E(j,i),VrhoA(j,i), &
                           VrhoB(j,i),VsigmaA(j,i),VsigmaB(j,i),VsigmaAB(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
    if( method == "b3pw91" .or. method == "B3PW91" .or. method == "rob3pw91" ) then
       do i = 1,M
          do j = 1,end_pot(i)

             x = gridx(j,i) ; y = gridy(j,i) ; z = gridz(j,i)
             call Cal_Rho_grid(x,y,z,rhoA,Rho_gridA(j,i,:),DensityA,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             call Cal_Rho_grid(x,y,z,rhoB,Rho_gridB(j,i,:),DensityB,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,"Y")
             sigmaA = Rho_gridA(j,i,1)**2 + Rho_gridA(j,i,2)**2 + Rho_gridA(j,i,3)**2
             sigmaB = Rho_gridB(j,i,1)**2 + Rho_gridB(j,i,2)**2 + Rho_gridB(j,i,3)**2
             sigmaAB = Rho_gridB(j,i,1)*Rho_gridA(j,i,1) + Rho_gridB(j,i,2)*Rho_gridA(j,i,2) + &
                       Rho_gridB(j,i,3)*Rho_gridA(j,i,3)
             call b3pw91_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E(j,i),VrhoA(j,i), &
                           VrhoB(j,i),VsigmaA(j,i),VsigmaB(j,i),VsigmaAB(j,i))
              
          end do
       end do
    end if
    !----------------------------------------------------------------
         
end subroutine rho_Func_Open
!=======================================================================
subroutine Cal_D_Phi(x,y,z,a,Ra,La,ATemp0,A3Temp0,YN)
    ! calculate the value of Phi and ▽Phi on every list.
    implicit none
    integer , parameter :: Node = 8760
    integer La(3)
    real(kind=8) r1a,expa,x2,y2,z2
    real(kind=8) x,y,z,Temp1,Temp,BNa
    real(kind=8) a,Ra(3)
    real(kind=8) ATemp0,A3Temp0(3)
    character(len=1) YN
    
    ATemp0 = 0.0d0 ; A3Temp0 = 0.0d0
    call Bas_Normal(a,La(1),La(2),La(3),BNa)
        
    x2 = (x-Ra(1))**La(1)
    y2 = (y-Ra(2))**La(2)
    z2 = (z-Ra(3))**La(3)
    r1a = (x-Ra(1))**2 + (y-Ra(2))**2 + (z-Ra(3))**2 
    expa = dexp(-a*r1a)
    ATemp0 = BNa*x2*y2*z2*expa
    
    if( YN == "N" ) return
    
    Temp1 = -2.0d0*a*ATemp0 
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( La(1) == 0 ) Temp = 0.0d0
    if ( La(1) >= 1 ) Temp = La(1)*BNa*(x-Ra(1))**(La(1)-1)*y2*z2*expa
    A3Temp0(1) = (x-Ra(1))*Temp1 + Temp
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( La(2) == 0 ) Temp = 0.0d0
    if ( La(2) >= 1 ) Temp = La(2)*BNa*x2*(y-Ra(2))**(La(2)-1)*z2*expa
    A3Temp0(2) = (y-Ra(2))*Temp1 + Temp
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if ( La(3) == 0 ) Temp = 0.0d0
    if ( La(3) >= 1 ) Temp = La(3)*BNa*x2*y2*(z-Ra(3))**(La(3)-1)*expa
    A3Temp0(3) = (z-Ra(3))*Temp1 + Temp

end subroutine Cal_D_Phi
!=======================================================================
subroutine Int_Func_E_XC(Node,M,Func_E,W_i,end_pot,Int_Res)
    ! exchange-correlation energy.
    implicit none
    integer Node
    integer M,i,j
    integer end_pot(M)
    real(kind=8) Int_Res,Func_E(Node,M)
    real(kind=8) W_i(Node,M)
    
    Int_Res = 0.0d0
    do i = 1,M

        do j = 1,end_pot(i)
                  
           Int_Res = Int_Res + W_i(j,i)*Func_E(j,i)

        end do
        
    end do  
    
end subroutine Int_Func_E_XC
!=======================================================================
