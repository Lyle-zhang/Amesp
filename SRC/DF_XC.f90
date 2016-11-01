subroutine Cal_Rho_grid(x,y,z,rho,Rho_grid,Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num,N,cont)
    ! Calculate Rho and ▽Rho.
    ! cont = "N" : Calculate Rho only ;
    ! cont = "Y" : Calculate Rho and ▽Rho ;
    implicit none
    integer N,ii,k,l,bas_num(N),NAtom(N),l_Basis(3,N)
    real(kind=8) x,y,z,rho
    real(kind=8) r2(N),x2(N),y2(N),z2(N),CBE(6,N)
    real(kind=8) Xi(N),DXi_x(N),DXi_y(N),DXi_z(N),Rho_grid(3)
    real(kind=8) Density(N,N),c_Basis(6,N,20),a_Basis(6,N,20),Rt(3,N) 
    real(kind=8) Temp1,BN,R_Tempx,R_Tempy,R_Tempz,Templ,Tempr
    character(len=1) cont
    
    !-------------------------------------------------------------------
    do ii = 1,N
       r2(ii) = (x-Rt(1,ii))**2+(y-Rt(2,ii))**2+(z-Rt(3,ii))**2
       x2(ii) = (x-Rt(1,ii))**l_Basis(1,ii)
       y2(ii) = (y-Rt(2,ii))**l_Basis(2,ii)
       z2(ii) = (z-Rt(3,ii))**l_Basis(3,ii) 
       
       Xi(ii) = 0.0d0
       do k = 1,bas_num(ii)
           call Bas_Normal(a_Basis(k,ii,NAtom(ii)),l_Basis(1,ii),l_Basis(2,ii),l_Basis(3,ii),BN)
           CBE(k,ii) = c_Basis(k,ii,NAtom(ii))*BN*dexp(-a_Basis(k,ii,NAtom(ii))*r2(ii))
           Xi(ii) = CBE(k,ii) + Xi(ii)
       end do
       Xi(ii) = Xi(ii)*x2(ii)*y2(ii)*z2(ii)
    end do 

    rho = 0.0d0         
    do k = 1,N
       Temp1 = 0.0d0
       do l = 1,N
          Temp1 = Temp1 + Density(k,l)*Xi(k)*Xi(l)
       end do
       rho = rho + Temp1
    end do
    
    if (cont == "N") return
    !-------------------------------------------------------------------
    do ii = 1,N
       DXi_x(ii) = 0.0d0
       DXi_y(ii) = 0.0d0
       DXi_z(ii) = 0.0d0

       Templ = x2(ii)*y2(ii)*z2(ii)
       do k = 1,bas_num(ii)
           
           Tempr = -2.0d0*a_Basis(k,ii,NAtom(ii))*Templ*CBE(k,ii)
           !------------------------------------------------------------
           if ( l_Basis(1,ii) == 0 ) then
               R_Tempx = 0.0d0
           else
               R_Tempx = CBE(k,ii)*l_Basis(1,ii)*((x-Rt(1,ii)))**(l_Basis(1,ii)-1) &
                    *y2(ii)*z2(ii)
           end if
           
           DXi_x(ii) = Tempr*((x-Rt(1,ii))) + R_Tempx + DXi_x(ii)
           !------------------------------------------------------------        
           if ( l_Basis(2,ii) == 0 ) then
               R_Tempy = 0.0d0
           else
               R_Tempy = CBE(k,ii)*l_Basis(2,ii)*x2(ii) &
                    *(y-Rt(2,ii))**(l_Basis(2,ii)-1)*z2(ii)
           end if
  
           DXi_y(ii) = Tempr*(y-Rt(2,ii)) + R_Tempy + DXi_y(ii)
           !------------------------------------------------------------
           if ( l_Basis(3,ii) == 0 ) then
               R_Tempz = 0.0d0
           else
               R_Tempz = CBE(k,ii)*l_Basis(3,ii)*x2(ii) &
                    *y2(ii)*(z-Rt(3,ii))**(l_Basis(3,ii)-1)
           end if
  
           DXi_z(ii) = Tempr*(z-Rt(3,ii)) + R_Tempz + DXi_z(ii)
                    
       end do
    end do    
    !-------------------------------------------------------------------
    Rho_grid(1) = 0.0d0         
    do k = 1,N
       Temp1 = 0.0d0
       do l = 1,N
          Temp1 = Temp1 + Density(k,l)*(DXi_x(k)*Xi(l) + Xi(k)*DXi_x(l))
       end do
       Rho_grid(1) = Rho_grid(1) + Temp1
    end do
    
    Rho_grid(2) = 0.0d0         
    do k = 1,N
       Temp1 = 0.0d0
       do l = 1,N
          Temp1 = Temp1 + Density(k,l)*(DXi_y(k)*Xi(l) + Xi(k)*DXi_y(l))
       end do
       Rho_grid(2) = Rho_grid(2) + Temp1
    end do
    
    Rho_grid(3) = 0.0d0         
    do k = 1,N
       Temp1 = 0.0d0
       do l = 1,N
          Temp1 = Temp1 + Density(k,l)*(DXi_z(k)*Xi(l) + Xi(k)*DXi_z(l))
       end do
       Rho_grid(3) = Rho_grid(3) + Temp1
    end do
    
end subroutine Cal_Rho_grid
!=======================================================================
!                                                                      !
!                 Vrho = delta{Func_E}/delta{rho}                      !
!                                                                      !
!=======================================================================
subroutine XAlpha_Close(rho,Func_E,Vrho)

    implicit none
    real(kind=8),parameter:: alpha =  0.70d0
    real(kind=8),parameter:: Coeff = -1.8610514726982d0
    real(kind=8),parameter:: Power =  0.33333333333333d0
    real(kind=8) rho,Vrho,Func_E
    
    Vrho = alpha*Coeff*(rho)**(Power) 
    Func_E = 2.0d0*0.75d0*Vrho*rho
    
end subroutine XAlpha_Close
!=======================================================================
subroutine XAlpha_Open(rhoA,rhoB,Func_E,VrhoA,VrhoB)

    implicit none
    real(kind=8),parameter:: alpha =  0.70d0
    real(kind=8),parameter:: Coeff = -1.8610514726982d0
    real(kind=8),parameter:: Power1 = 1.33333333333333d0
    real(kind=8),parameter:: Power2 = 0.33333333333333d0
    real(kind=8) rhoA,rhoB,Func_E,VrhoA,VrhoB
    
    VrhoA = alpha*Coeff*((rhoA)**(Power2))
    VrhoB = alpha*Coeff*((rhoB)**(Power2))
    Func_E = 0.75d0*VrhoA*rhoA + 0.75d0*VrhoB*rhoB
    
end subroutine XAlpha_Open
!=======================================================================
subroutine LSDA_Close(rho,Func_E,Vrho)

    implicit none

    real(kind=8) rho,Vrho,Func_E
    real(kind=8) rhoa1(1),sigmaaa1(1),zk_x(1),vrhoa_x(1),zk_c(1),vrhoa_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmaaa_c(1),v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0
    call rks_x_lda(1,1,rhoa1,sigmaaa1,zk_x,vrhoa_x,vsigmaaa_x,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_vwn5(1,1,rhoa1,sigmaaa1,zk_c,vrhoa_c,vsigmaaa_c,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    Func_E = zk_x(1) + zk_c(1)
    Vrho = vrhoa_x(1) + vrhoa_c(1)
    
end subroutine LSDA_Close
!=======================================================================
subroutine LSDA_Open(rhoA,rhoB,Func_E,VrhoA,VrhoB)

    implicit none
    real(kind=8) rhoA,rhoB,Func_E,VrhoA,VrhoB
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1),zk_x(1),zk_c(1)
    real(kind=8) vrhoa_x(1),vrhob_x(1),vrhoa_c(1),vrhob_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmabb_x(1),vsigmaab_x(1),vsigmaaa_c(1),vsigmabb_c(1),vsigmaab_c(1)
    real(kind=8) v2rhoa2(1),v2rhob2(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    
    call uks_x_lda(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_x, &
         vrhoa_x,vrhob_x,vsigmaaa_x,vsigmabb_x,vsigmaab_x,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)

    call uks_c_vwn5(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_c, &
         vrhoa_c,vrhob_c,vsigmaaa_c,vsigmabb_c,vsigmaab_c,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
                 
    Func_E = zk_x(1) + zk_c(1)
    VrhoA = vrhoa_x(1) + vrhoa_c(1) 
    VrhoB = vrhob_x(1) + vrhob_c(1)
         
end subroutine LSDA_Open
!=======================================================================
subroutine pbe_Close(rho,sigma,Func_E,Vrho,Vsigma)

    implicit none

    real(kind=8) rho,sigma,Vrho,Vsigma,Func_E
    real(kind=8) rhoa1(1),sigmaaa1(1),zk_x(1),vrhoa_x(1),zk_c(1),vrhoa_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmaaa_c(1),v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0 ; sigmaaa1 = sigma*4.0d0
    call rks_x_pbe(1,1,rhoa1,sigmaaa1,zk_x,vrhoa_x,vsigmaaa_x,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_pbe(1,1,rhoa1,sigmaaa1,zk_c,vrhoa_c,vsigmaaa_c,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    Func_E = zk_x(1) + zk_c(1)
    Vrho = vrhoa_x(1) + vrhoa_c(1)
    Vsigma = vsigmaaa_x(1) + vsigmaaa_c(1)
    Vsigma = Vsigma/2.0d0
    
end subroutine pbe_Close
!=======================================================================
subroutine pbe_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)

    implicit none
    real(kind=8) rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1),zk_x(1),zk_c(1)
    real(kind=8) vrhoa_x(1),vrhob_x(1),vrhoa_c(1),vrhob_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmabb_x(1),vsigmaab_x(1),vsigmaaa_c(1),vsigmabb_c(1),vsigmaab_c(1)
    real(kind=8) v2rhoa2(1),v2rhob2(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    sigmaaa1 = sigmaA ; sigmabb1 = sigmaB ; sigmaab1 = sigmaAB
    
    call uks_x_pbe(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_x, &
         vrhoa_x,vrhob_x,vsigmaaa_x,vsigmabb_x,vsigmaab_x,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)

    call uks_c_pbe(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_c, &
         vrhoa_c,vrhob_c,vsigmaaa_c,vsigmabb_c,vsigmaab_c,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
                 
    Func_E = zk_x(1) + zk_c(1)
    VrhoA = vrhoa_x(1) + vrhoa_c(1) 
    VrhoB = vrhob_x(1) + vrhob_c(1)
    VsigmaA = vsigmaaa_x(1) + vsigmaaa_c(1)
    VsigmaB = vsigmabb_x(1) + vsigmabb_c(1)
    VsigmaAB = vsigmaab_x(1) + vsigmaab_c(1) 
         
end subroutine pbe_Open
!=======================================================================
subroutine pw91_Close(rho,sigma,Func_E,Vrho,Vsigma)

    implicit none

    real(kind=8) rho,sigma,Vrho,Vsigma,Func_E
    real(kind=8) rhoa1(1),sigmaaa1(1),zk_x(1),vrhoa_x(1),zk_c(1),vrhoa_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmaaa_c(1),v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0 ; sigmaaa1 = sigma*4.0d0
    call rks_x_pw91(1,1,rhoa1,sigmaaa1,zk_x,vrhoa_x,vsigmaaa_x,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_pw91(1,1,rhoa1,sigmaaa1,zk_c,vrhoa_c,vsigmaaa_c,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    Func_E = zk_x(1) + zk_c(1)
    Vrho = vrhoa_x(1) + vrhoa_c(1)
    Vsigma = vsigmaaa_x(1) + vsigmaaa_c(1)
    Vsigma = Vsigma/2.0d0
    
end subroutine pw91_Close
!=======================================================================
subroutine pw91_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)

    implicit none
    real(kind=8) rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1),zk_x(1),zk_c(1)
    real(kind=8) vrhoa_x(1),vrhob_x(1),vrhoa_c(1),vrhob_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmabb_x(1),vsigmaab_x(1),vsigmaaa_c(1),vsigmabb_c(1),vsigmaab_c(1)
    real(kind=8) v2rhoa2(1),v2rhob2(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    sigmaaa1 = sigmaA ; sigmabb1 = sigmaB ; sigmaab1 = sigmaAB
    
    call uks_x_pw91(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_x, &
         vrhoa_x,vrhob_x,vsigmaaa_x,vsigmabb_x,vsigmaab_x,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)

    call uks_c_pw91(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_c, &
         vrhoa_c,vrhob_c,vsigmaaa_c,vsigmabb_c,vsigmaab_c,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
                 
    Func_E = zk_x(1) + zk_c(1)
    VrhoA = vrhoa_x(1) + vrhoa_c(1) 
    VrhoB = vrhob_x(1) + vrhob_c(1)
    VsigmaA = vsigmaaa_x(1) + vsigmaaa_c(1)
    VsigmaB = vsigmabb_x(1) + vsigmabb_c(1)
    VsigmaAB = vsigmaab_x(1) + vsigmaab_c(1) 
         
end subroutine pw91_Open
!=======================================================================
subroutine B3LYP_Close(rho,sigma,Func_E,Vrho,Vsigma)

    implicit none

    real(kind=8) rho,sigma,Vrho,Vsigma,Func_E
    real(kind=8) a0,ax,ac
    real(kind=8) lda_x(1),b88_x(1),vrhoa_lda_x(1),vrhoa_b88_x(1)
    real(kind=8) vsigmaaa_lda_x(1),vsigmaaa_b88_x(1),vsigmaaa_vwn5_c(1),vsigmaaa_lyp_c(1)
    real(kind=8) vwn5_c(1),lyp_c(1),vrhoa_vwn5_c(1),vrhoa_lyp_c(1)
    real(kind=8) rhoa1(1),sigmaaa1(1)
    real(kind=8) v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0 ; sigmaaa1 = sigma*4.0d0
    call rks_x_lda(1,1,rhoa1,sigmaaa1,lda_x,vrhoa_lda_x,vsigmaaa_lda_x, &
                   v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_x_b88(1,1,rhoa1,sigmaaa1,b88_x,vrhoa_b88_x,vsigmaaa_b88_x, &
                   v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_vwn5(1,1,rhoa1,sigmaaa1,vwn5_c,vrhoa_vwn5_c,vsigmaaa_vwn5_c, &
                    v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_lyp(1,1,rhoa1,sigmaaa1,lyp_c,vrhoa_lyp_c,vsigmaaa_lyp_c, &
                   v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    a0 = 0.2d0 ; ax = 0.72d0 ; ac = 0.81d0
    
    Func_E = (1.0d0-a0-ax)*lda_x(1) + ax*b88_x(1) + (1.0d0-ac)*vwn5_c(1) + ac*lyp_c(1)
    Vrho = (1.0d0-a0-ax)*vrhoa_lda_x(1) + ax*vrhoa_b88_x(1) + &
           (1.0d0-ac)*vrhoa_vwn5_c(1) + ac*vrhoa_lyp_c(1)
    Vsigma = ax*vsigmaaa_b88_x(1) + &
           (1.0d0-ac)*vsigmaaa_vwn5_c(1) + ac*vsigmaaa_lyp_c(1)
    Vsigma = Vsigma/2.0d0
   
end subroutine B3LYP_Close
!=======================================================================
subroutine B3LYP_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)

    implicit none

    real(kind=8) rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB
    real(kind=8) a0,ax,ac
    real(kind=8) lda_x(1),b88_x(1),vwn5_c(1),lyp_c(1)
    real(kind=8) vrhoa_lda_x(1),vrhob_lda_x(1),vsigmaaa_lda_x(1),vsigmabb_lda_x(1)
    real(kind=8) vrhoa_b88_x(1),vrhob_b88_x(1),vsigmaaa_b88_x(1),vsigmabb_b88_x(1)
    real(kind=8) vrhoa_vwn5_c(1),vrhob_vwn5_c(1),vsigmaaa_vwn5_c(1),vsigmabb_vwn5_c(1)
    real(kind=8) vrhoa_lyp_c(1),vrhob_lyp_c(1),vsigmaaa_lyp_c(1),vsigmabb_lyp_c(1) 
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1)
    real(kind=8) vsigmaab_lda_x(1),vsigmaab_b88_x(1),v2rhoa2(1),v2rhob2(1)
    real(kind=8) vsigmaab_vwn5_c(1),vsigmaab_lyp_c(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    sigmaaa1 = sigmaA ; sigmabb1 = sigmaB ; sigmaab1 = sigmaAB
    call uks_x_lda(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,lda_x, &
                   vrhoa_lda_x,vrhob_lda_x,vsigmaaa_lda_x,vsigmabb_lda_x &
                   ,vsigmaab_lda_x,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    call uks_x_b88(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,b88_x, &
                   vrhoa_b88_x,vrhob_b88_x,vsigmaaa_b88_x,vsigmabb_b88_x &
                   ,vsigmaab_b88_x,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    call uks_c_vwn5(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,vwn5_c, &
                   vrhoa_vwn5_c,vrhob_vwn5_c,vsigmaaa_vwn5_c,vsigmabb_vwn5_c &
                   ,vsigmaab_vwn5_c,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    call uks_c_lyp(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,lyp_c, &
                   vrhoa_lyp_c,vrhob_lyp_c,vsigmaaa_lyp_c,vsigmabb_lyp_c &
                   ,vsigmaab_lyp_c,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    
    a0 = 0.2d0 ; ax = 0.72d0 ; ac = 0.81d0
    
    Func_E = (1.0d0-a0-ax)*lda_x(1) + ax*b88_x(1) + (1.0d0-ac)*vwn5_c(1) + ac*lyp_c(1) 
    VrhoA = (1.0d0-a0-ax)*vrhoa_lda_x(1) + ax*vrhoa_b88_x(1) + &
            (1.0d0-ac)*vrhoa_vwn5_c(1) + ac*vrhoa_lyp_c(1)
    VrhoB = (1.0d0-a0-ax)*vrhob_lda_x(1) + ax*vrhob_b88_x(1) + &
            (1.0d0-ac)*vrhob_vwn5_c(1) + ac*vrhob_lyp_c(1)
    VsigmaA = ax*vsigmaaa_b88_x(1) + &
              (1.0d0-ac)*vsigmaaa_vwn5_c(1) + ac*vsigmaaa_lyp_c(1)
    VsigmaB = ax*vsigmabb_b88_x(1) + &
              (1.0d0-ac)*vsigmabb_vwn5_c(1) + ac*vsigmabb_lyp_c(1)
    VsigmaAB = ax*vsigmaab_b88_x(1) + &
              (1.0d0-ac)*vsigmaab_vwn5_c(1) + ac*vsigmaab_lyp_c(1)
   
end subroutine B3LYP_Open
!=======================================================================
subroutine blyp_Close(rho,sigma,Func_E,Vrho,Vsigma)

    implicit none

    real(kind=8) rho,sigma,Vrho,Vsigma,Func_E
    real(kind=8) rhoa1(1),sigmaaa1(1),zk_x(1),vrhoa_x(1),zk_c(1),vrhoa_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmaaa_c(1),v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0 ; sigmaaa1 = sigma*4.0d0
    call rks_x_b88(1,1,rhoa1,sigmaaa1,zk_x,vrhoa_x,vsigmaaa_x,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_lyp(1,1,rhoa1,sigmaaa1,zk_c,vrhoa_c,vsigmaaa_c,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    Func_E = zk_x(1) + zk_c(1)
    Vrho = vrhoa_x(1) + vrhoa_c(1)
    Vsigma = vsigmaaa_x(1) + vsigmaaa_c(1)
    Vsigma = Vsigma/2.0d0
    
end subroutine blyp_Close
!=======================================================================
subroutine blyp_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)

    implicit none
    real(kind=8) rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1),zk_x(1),zk_c(1)
    real(kind=8) vrhoa_x(1),vrhob_x(1),vrhoa_c(1),vrhob_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmabb_x(1),vsigmaab_x(1),vsigmaaa_c(1),vsigmabb_c(1),vsigmaab_c(1)
    real(kind=8) v2rhoa2(1),v2rhob2(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    sigmaaa1 = sigmaA ; sigmabb1 = sigmaB ; sigmaab1 = sigmaAB
    
    call uks_x_b88(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_x, &
         vrhoa_x,vrhob_x,vsigmaaa_x,vsigmabb_x,vsigmaab_x,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)

    call uks_c_lyp(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_c, &
         vrhoa_c,vrhob_c,vsigmaaa_c,vsigmabb_c,vsigmaab_c,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
                 
    Func_E = zk_x(1) + zk_c(1)
    VrhoA = vrhoa_x(1) + vrhoa_c(1) 
    VrhoB = vrhob_x(1) + vrhob_c(1)
    VsigmaA = vsigmaaa_x(1) + vsigmaaa_c(1)
    VsigmaB = vsigmabb_x(1) + vsigmabb_c(1)
    VsigmaAB = vsigmaab_x(1) + vsigmaab_c(1) 
         
end subroutine blyp_Open
!=======================================================================
subroutine pbe0_Close(rho,sigma,Func_E,Vrho,Vsigma)

    implicit none

    real(kind=8) rho,sigma,Vrho,Vsigma,Func_E
    real(kind=8) rhoa1(1),sigmaaa1(1),zk_x(1),vrhoa_x(1),zk_c(1),vrhoa_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmaaa_c(1),v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0 ; sigmaaa1 = sigma*4.0d0
    call rks_x_pbe(1,1,rhoa1,sigmaaa1,zk_x,vrhoa_x,vsigmaaa_x,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_pbe(1,1,rhoa1,sigmaaa1,zk_c,vrhoa_c,vsigmaaa_c,v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    Func_E = zk_x(1)*0.75d0 + zk_c(1)
    Vrho = vrhoa_x(1)*0.75d0 + vrhoa_c(1)
    Vsigma = vsigmaaa_x(1)*0.75d0 + vsigmaaa_c(1)
    Vsigma = Vsigma/4.0d0
    
end subroutine pbe0_Close
!=======================================================================
subroutine pbe0_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)

    implicit none
    real(kind=8) rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1),zk_x(1),zk_c(1)
    real(kind=8) vrhoa_x(1),vrhob_x(1),vrhoa_c(1),vrhob_c(1)
    real(kind=8) vsigmaaa_x(1),vsigmabb_x(1),vsigmaab_x(1),vsigmaaa_c(1),vsigmabb_c(1),vsigmaab_c(1)
    real(kind=8) v2rhoa2(1),v2rhob2(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    sigmaaa1 = sigmaA ; sigmabb1 = sigmaB ; sigmaab1 = sigmaAB
    
    call uks_x_pbe(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_x, &
         vrhoa_x,vrhob_x,vsigmaaa_x,vsigmabb_x,vsigmaab_x,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)

    call uks_c_pbe(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,zk_c, &
         vrhoa_c,vrhob_c,vsigmaaa_c,vsigmabb_c,vsigmaab_c,v2rhoa2,v2rhob2, &
         v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
         v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
         v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
                 
    Func_E = zk_x(1)*0.75d0 + zk_c(1)
    VrhoA = vrhoa_x(1)*0.75d0 + vrhoa_c(1) 
    VrhoB = vrhob_x(1)*0.75d0 + vrhob_c(1)
    VsigmaA = vsigmaaa_x(1)*0.75d0 + vsigmaaa_c(1)
    VsigmaB = vsigmabb_x(1)*0.75d0 + vsigmabb_c(1)
    VsigmaAB = vsigmaab_x(1)*0.75d0 + vsigmaab_c(1) 
         
end subroutine pbe0_Open
!=======================================================================
subroutine B3PW91_Close(rho,sigma,Func_E,Vrho,Vsigma)

    implicit none

    real(kind=8) rho,sigma,Vrho,Vsigma,Func_E
    real(kind=8) a0,ax,ac
    real(kind=8) lda_x(1),b88_x(1),vrhoa_lda_x(1),vrhoa_b88_x(1)
    real(kind=8) vsigmaaa_lda_x(1),vsigmaaa_b88_x(1),vsigmaaa_vwn5_c(1),vsigmaaa_pw91_c(1)
    real(kind=8) vwn5_c(1),pw91_c(1),vrhoa_vwn5_c(1),vrhoa_pw91_c(1)
    real(kind=8) rhoa1(1),sigmaaa1(1)
    real(kind=8) v2rhoa2(1),v2rhoasigmaaa(1),v2sigmaaa2(1)
    
    rhoa1 = rho*2.0d0 ; sigmaaa1 = sigma*4.0d0
    call rks_x_lda(1,1,rhoa1,sigmaaa1,lda_x,vrhoa_lda_x,vsigmaaa_lda_x, &
                   v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_x_b88(1,1,rhoa1,sigmaaa1,b88_x,vrhoa_b88_x,vsigmaaa_b88_x, &
                   v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_vwn5(1,1,rhoa1,sigmaaa1,vwn5_c,vrhoa_vwn5_c,vsigmaaa_vwn5_c, &
                    v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    call rks_c_pw91(1,1,rhoa1,sigmaaa1,pw91_c,vrhoa_pw91_c,vsigmaaa_pw91_c, &
                   v2rhoa2,v2rhoasigmaaa,v2sigmaaa2)
    
    a0 = 0.2d0 ; ax = 0.72d0 ; ac = 0.81d0
    
    Func_E = (1.0d0-a0-ax)*lda_x(1) + ax*b88_x(1) + (1.0d0-ac)*vwn5_c(1) + ac*pw91_c(1)
    Vrho = (1.0d0-a0-ax)*vrhoa_lda_x(1) + ax*vrhoa_b88_x(1) + &
           (1.0d0-ac)*vrhoa_vwn5_c(1) + ac*vrhoa_pw91_c(1)
    Vsigma = ax*vsigmaaa_b88_x(1) + &
           (1.0d0-ac)*vsigmaaa_vwn5_c(1) + ac*vsigmaaa_pw91_c(1)
    Vsigma = Vsigma/2.0d0
   
end subroutine B3PW91_Close
!=======================================================================
subroutine B3PW91_Open(rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB)

    implicit none

    real(kind=8) rhoA,rhoB,sigmaA,sigmaB,sigmaAB,Func_E,VrhoA,VrhoB,VsigmaA,VsigmaB,VsigmaAB
    real(kind=8) a0,ax,ac
    real(kind=8) lda_x(1),b88_x(1),vwn5_c(1),pw91_c(1)
    real(kind=8) vrhoa_lda_x(1),vrhob_lda_x(1),vsigmaaa_lda_x(1),vsigmabb_lda_x(1)
    real(kind=8) vrhoa_b88_x(1),vrhob_b88_x(1),vsigmaaa_b88_x(1),vsigmabb_b88_x(1)
    real(kind=8) vrhoa_vwn5_c(1),vrhob_vwn5_c(1),vsigmaaa_vwn5_c(1),vsigmabb_vwn5_c(1)
    real(kind=8) vrhoa_pw91_c(1),vrhob_pw91_c(1),vsigmaaa_pw91_c(1),vsigmabb_pw91_c(1) 
    real(kind=8) rhoa1(1),rhob1(1),sigmaaa1(1),sigmabb1(1),sigmaab1(1)
    real(kind=8) vsigmaab_lda_x(1),vsigmaab_b88_x(1),v2rhoa2(1),v2rhob2(1)
    real(kind=8) vsigmaab_vwn5_c(1),vsigmaab_pw91_c(1)
    real(kind=8) v2rhoab(1),v2rhoasigmaaa(1),v2rhoasigmaab(1),v2rhoasigmabb(1)
    real(kind=8) v2rhobsigmabb(1),v2rhobsigmaab(1),v2rhobsigmaaa(1),v2sigmaaa2(1)
    real(kind=8) v2sigmaaaab(1),v2sigmaaabb(1),v2sigmaab2(1),v2sigmaabbb(1),v2sigmabb2(1)
    
    rhoa1 = rhoA ; rhob1 = rhoB
    sigmaaa1 = sigmaA ; sigmabb1 = sigmaB ; sigmaab1 = sigmaAB
    call uks_x_lda(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,lda_x, &
                   vrhoa_lda_x,vrhob_lda_x,vsigmaaa_lda_x,vsigmabb_lda_x &
                   ,vsigmaab_lda_x,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    call uks_x_b88(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,b88_x, &
                   vrhoa_b88_x,vrhob_b88_x,vsigmaaa_b88_x,vsigmabb_b88_x &
                   ,vsigmaab_b88_x,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    call uks_c_vwn5(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,vwn5_c, &
                   vrhoa_vwn5_c,vrhob_vwn5_c,vsigmaaa_vwn5_c,vsigmabb_vwn5_c &
                   ,vsigmaab_vwn5_c,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    call uks_c_pw91(1,1,rhoa1,rhob1,sigmaaa1,sigmabb1,sigmaab1,pw91_c, &
                   vrhoa_pw91_c,vrhob_pw91_c,vsigmaaa_pw91_c,vsigmabb_pw91_c &
                   ,vsigmaab_pw91_c,v2rhoa2,v2rhob2, &
                   v2rhoab,v2rhoasigmaaa,v2rhoasigmaab,v2rhoasigmabb, &
                   v2rhobsigmabb,v2rhobsigmaab,v2rhobsigmaaa,v2sigmaaa2, &
                   v2sigmaaaab,v2sigmaaabb,v2sigmaab2,v2sigmaabbb,v2sigmabb2)
    
    a0 = 0.2d0 ; ax = 0.72d0 ; ac = 0.81d0
    
    Func_E = (1.0d0-a0-ax)*lda_x(1) + ax*b88_x(1) + (1.0d0-ac)*vwn5_c(1) + ac*pw91_c(1) 
    VrhoA = (1.0d0-a0-ax)*vrhoa_lda_x(1) + ax*vrhoa_b88_x(1) + &
            (1.0d0-ac)*vrhoa_vwn5_c(1) + ac*vrhoa_pw91_c(1)
    VrhoB = (1.0d0-a0-ax)*vrhob_lda_x(1) + ax*vrhob_b88_x(1) + &
            (1.0d0-ac)*vrhob_vwn5_c(1) + ac*vrhob_pw91_c(1)
    VsigmaA = ax*vsigmaaa_b88_x(1) + &
              (1.0d0-ac)*vsigmaaa_vwn5_c(1) + ac*vsigmaaa_pw91_c(1)
    VsigmaB = ax*vsigmabb_b88_x(1) + &
              (1.0d0-ac)*vsigmabb_vwn5_c(1) + ac*vsigmabb_pw91_c(1)
    VsigmaAB = ax*vsigmaab_b88_x(1) + &
              (1.0d0-ac)*vsigmaab_vwn5_c(1) + ac*vsigmaab_pw91_c(1)
   
end subroutine B3PW91_Open
!=======================================================================
