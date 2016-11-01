!=======================================================================
!                             ERI (sp)                                 !
!     See : Explicit formulae for integrals of s and p type GTFs       ! 
!           ---- CAROL A. BAXTER AND DAVID B. COOK *                   ! 
!=======================================================================
subroutine Two_Electron_Int(Two_Electron,EI,N,Nmem,Neri,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,t2,d2)

    implicit none
    integer N,Nmem,Neri,i,j,k,l,bas_num(N),NAtom(N)
    integer EI(N,N),l_Basis(3,N)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2),Rt(3,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) Temp1
    character(len=10) t2
    character(len=8) d2
    
    k = 0
    do i = 1,N
       do j = 1,i
          k = K + 1
          EI(j,i) = k
          EI(i,j) = k
       end do
    end do 

    
    Two_Electron = 0.0d0
    !-------------------------Schwars Screening-------------------------
    if ( N > Nmem ) then
       do i = 1,N
          do j = i,N  
             call ERI(N,i,j,i,j,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,Temp1)
             Two_Electron(EI(j,i),1) = Temp1   !Two_Electron
          end do
       end do
    else
       do i = 1,N
          do j = i,N  
             call ERI(N,i,j,i,j,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,Temp1)
             Two_Electron(EI(j,i),EI(j,i)) = Temp1   !Two_Electron
          end do
       end do
    end if
    !-------------------------------------------------------------------
    if ( N > Nmem ) then
       call date_and_time( date = d2 , time = t2)
       open(unit=1110,file="/tmp/"//d2//t2//".amesp",form="unformatted")
        
       do i = 1,N
         do j = 1,i
           do k = 1,N
             do l = 1,k
       
                if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                if(dsqrt(dabs(Two_Electron(EI(j,i),1)*Two_Electron(EI(l,k),1))) < 1.0d-8 ) cycle
        
                call ERI(N,i,j,k,l,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,Temp1)
                write(1110) Temp1
                  
             end do
           end do
         end do
       end do    
       close(1110)  
    else !--------------------------------------------------------------
       do i = 1,N
         do j = i,N
           do k = 1,N
             do l = k,N
       
               if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
               if( i == k .and. j == l ) cycle
               if(dsqrt(dabs(Two_Electron(EI(j,i),EI(j,i))*Two_Electron(EI(l,k),EI(l,k)))) < 1.0d-8 ) cycle
        
               call ERI(N,i,j,k,l,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,Temp1)
               Two_Electron(EI(l,k),EI(j,i)) = Temp1   !Two_Electron
               Two_Electron(EI(j,i),EI(l,k)) = Temp1
             end do
           end do
         end do
       end do
    end if
    !-------------------------------------------------------------------
end subroutine Two_Electron_Int
!=======================================================================
subroutine ERI(N,i,j,k,l,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom,Temp1)

    implicit none
    integer N,ii,i,j,k,l,i1,i2,i3,i4,bas_num(N),NAtom(N)
    integer l_Basis(3,N),Delta(3,3)
    real(kind=8) Rt(3,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) Temp0,Temp1,Temp2,Temp3,Temp4 

    
    Delta = 0
    do ii = 1,3
       Delta(ii,ii) = 1
    end do  
    
    Temp1 = 0.0d0
    do i1 = 1,bas_num(i)
      Temp2 = 0.0d0
      do i2 = 1,bas_num(j)
        Temp3 = 0.0d0
        do i3 = 1,bas_num(k)
          Temp4 = 0.0d0
          do i4 = 1,bas_num(l)
            !call E_E_Repulsion_sp(a,b,c,d,Ra,Rb,Rc,La,Lb,Lc,Ld,Ld,Delta,S)     
             call E_E_Repulsion_sp(a_Basis(i1,i,NAtom(i)),a_Basis(i2,j,NAtom(j)),a_Basis(i3,k,NAtom(k)),&
                  a_Basis(i4,l,NAtom(l)),Rt(:,i),Rt(:,j),Rt(:,k),Rt(:,l),l_Basis(:,i),l_Basis(:,j), &
                  l_Basis(:,k),l_Basis(:,l),Delta,Temp0)
                  Temp4 = Temp4 + Temp0*c_Basis(i1,i,NAtom(i))*c_Basis(i2,j,NAtom(j))* &
                          c_Basis(i3,k,NAtom(k))*c_Basis(i4,l,NAtom(l))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do  

end subroutine ERI
!=======================================================================
subroutine E_E_Repulsion_sp(a,b,c,d,Ra,Rb,Rc,Rd,La,Lb,Lc,Ld,Delta,S)

    implicit none
    integer i,j,k,l,La(3),Lb(3),Lc(3),Ld(3)
    integer ijkl,Delta(3,3),sa,sb,sc,sd
    real(kind=8),parameter::Pi = 3.141592653589793d0
    real(kind=8) a,b,c,d,AB,CD,S
    real(kind=8) K1,k2,Rp(3),Rq(3)
    real(kind=8) S1,S2,S4,LAMBDA
    real(kind=8) BNabcd
    real(kind=8) S00ab,S00cd,Si0ab,S0jab,Sk0cd,S0lcd,Sijab,Sklcd
    real(kind=8) G_i000,G_00k0,G_i0k0,G_0j00,G_ij00
    real(kind=8) G_0jk0,G_000l,G_i00l,G_0j0l,G_00kl
    real(kind=8) G_0jkl,G_i0kl,G_ijk0,G_ij0l,G_ijkl
    real(kind=8) PQ,w,F_0,F_1,F_2,F_3,F_4
    real(kind=8) Ra(3),Rb(3),Rc(3),Rd(3)
    
    S = 0.0d0 
    call GTO_s(a,b,Ra,Rb,K1,AB,Rp)
    call GTO_s(c,d,Rc,Rd,K2,CD,Rq)
    
    PQ = (Rp(1)-Rq(1))**2+(Rp(2)-Rq(2))**2+(Rp(3)-Rq(3))**2
    S1 = a + b ; S2 = c + d ; S4 = S1 + S2
    w = PQ*S1*S2/S4 ;    
    
    sa = sum(La(:)) ; sb = sum(Lb(:)) ; sc = sum(Lc(:)) ; sd = sum(Ld(:))
    ijkl = sa + sb + sc + sd
    !------------------------- <ss|1/rC|ss> ----------------------------
    if( ijkl == 0 ) then
        call F0(w,F_0)
        S = K1*K2*F_0*(9.0270333367641006d0*(a*b*c*d)**(0.75d0))/((S1)*(S2)* &
            dsqrt(S4))
        return
    end if
    !-------------------------------------------------------------------

    LAMBDA = 2*dsqrt(S1*S2/(S4*Pi))
      
    !------------------------- <ps|1/rC|ss> ----------------------------
    if( ijkl == 1 ) then
    
        if( sa == 1 ) then 
            !------<ps|1/rC|ss>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
            
            S = Si0ab*S00cd*F_0 + S00ab*S00cd*G_i000
            BNabcd = 0.51602455093119183d0*a**1.25*(b*c*d)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sb == 1 ) then 
            !------<sp|1/rC|ss>------!
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_i000)
            
            S = S0jab*S00cd*F_0 + S00ab*S00cd*G_i000
            BNabcd = 0.51602455093119183d0*b**1.25*(a*c*d)**0.75
            S = S*BNabcd*LAMBDA
            
            return          
        end if
        
        if( sc == 1 ) then 
            !------<ss|1/rC|ps>------!
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)
            
            S = S00ab*Sk0cd*F_0 + S00ab*S00cd*G_00k0
            BNabcd = 0.51602455093119183d0*c**1.25*(b*a*d)**0.75
            S = S*BNabcd*LAMBDA
            
            return           
        end if
        
        if( sd == 1 ) then 
            !------<ss|1/rC|sp>------!
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            
            S = S00ab*S0lcd*F_0 + S00ab*S00cd*G_000l
            BNabcd = 0.51602455093119183d0*d**1.25*(b*c*a)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
    end if
    !------------------------- <ps|1/rC|ps> ----------------------------
    if( ijkl == 2 ) then
    
        if( sa == 1 .and. sc == 1 ) then
            !------<ps|1/rC|ps>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)            
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
            call Gi0j0(i,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i0k0)
            
            S = Si0ab*S00cd*G_00k0 + Si0ab*Sk0cd*F_0 + S00ab*S00cd*G_i0k0 &
                + S00ab*Sk0cd*G_i000
            BNabcd = 1.03204910186238365d0*(a*c)**1.25*(b*d)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sb == 1 .and. sc == 1 ) then
            !------<sp|1/rC|ps>------!
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do
            
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)            
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00)
            call Gi0j0(j,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0jk0)
            
            S = S0jab*S00cd*G_00k0 + S0jab*Sk0cd*F_0 + S00ab*S00cd*G_0jk0 &
                + S00ab*Sk0cd*G_0j00
            BNabcd = 1.03204910186238365d0*(b*c)**1.25*(a*d)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sa == 1 .and. sd == 1 ) then
            !------<ps|1/rC|sp>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)            
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
            call Gi0j0(i,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i00l)
            
            S = Si0ab*S00cd*G_000l + Si0ab*S0lcd*F_0 + S00ab*S00cd*G_i00l &
                + S00ab*S0lcd*G_i000
            BNabcd = 1.03204910186238365d0*(a*d)**1.25*(b*c)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sb == 1 .and. sd == 1 ) then
            !------<sp|1/rC|sp>------!
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do
            
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)            
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00)
            call Gi0j0(j,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0j0l)
            
            S = S0jab*S00cd*G_000l + S0jab*S0lcd*F_0 + S00ab*S00cd*G_0j0l &
                + S00ab*S0lcd*G_0j00
            BNabcd = 1.03204910186238365d0*(b*d)**1.25*(a*c)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sa == 1 .and. sb == 1 ) then
            !------<pp|1/rC|ss>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)            
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            call Sabij(a,b,Ra,Rb,i,j,Delta,S00ab,Sijab)
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00)
            call Gij00(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_ij00)
            
            S = Sijab*S00cd*F_0 + Si0ab*S00cd*G_0j00 + S0jab*S00cd*G_i000 &
                + S00ab*S00cd*G_ij00
            BNabcd = 1.03204910186238365d0*(a*b)**1.25*(d*c)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sc == 1 .and. sd == 1 ) then
            !------<ss|1/rC|pp>------!
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do
            
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do
            
            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)            
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            call Sabij(c,d,Rc,Rd,k,l,Delta,S00cd,Sklcd)
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            call G00ij(k,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_00kl)
            
            S = Sklcd*S00ab*F_0 + Sk0cd*S00ab*G_000l + S0lcd*S00ab*G_00k0 &
                + S00ab*S00cd*G_00kl
            BNabcd = 1.03204910186238365d0*(c*d)**1.25*(a*b)**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
    end if
    !------------------------- <pp|1/rC|ps> ----------------------------
    if( ijkl == 3 ) then
    
        if( sd == 0 ) then 
            !------<pp|1/rC|ps>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do 
            
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do

            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2) ; call F3(w,F_3)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)  
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            call Sabij(a,b,Ra,Rb,i,j,Delta,S00ab,Sijab) 
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00)
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
            call Gij00(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_ij00)  
            
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)  
            call Gi0j0(j,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0jk0)  
            call Gi0j0(i,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i0k0)
            call Gijk0(i,j,k,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_ijk0)                      
           
            S = Sijab*( Sk0cd*F_0    + S00cd*G_00k0 ) + &
                Si0ab*( Sk0cd*G_0j00 + S00cd*G_0jk0 ) + &
                S0jab*( Sk0cd*G_i000 + S00cd*G_i0k0 ) + &
                S00ab*( Sk0cd*G_ij00 + S00cd*G_ijk0 ) 
                
            BNabcd = 2.06409820372476731d0*(a*b*c)**1.25*d**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
           
        if( sc == 0 ) then 
            !------<pp|1/rC|sp>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do 
            
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do

            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2) ; call F3(w,F_3)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)  
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            call Sabij(a,b,Ra,Rb,i,j,Delta,S00ab,Sijab) 
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00)
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
            call Gij00(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_ij00)  
            
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            call Gi0j0(j,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0j0l)
            call Gi0j0(i,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i00l)
            call Gijk0(i,j,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_ij0l)                      
           
            S = Sijab*( S0lcd*F_0    + S00cd*G_000l ) + &
                Si0ab*( S0lcd*G_0j00 + S00cd*G_0j0l ) + &
                S0jab*( S0lcd*G_i000 + S00cd*G_i00l ) + &
                S00ab*( S0lcd*G_ij00 + S00cd*G_ij0l ) 
                
            BNabcd = 2.06409820372476731d0*(a*b*d)**1.25*c**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if

        if( sb == 0 ) then 
            !------<ps|1/rC|pp>------!
            do i = 1,3
               if( La(i) == 1 ) exit
            end do
            
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do
            
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do

            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2) ; call F3(w,F_3)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)  
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            call Sabij(c,d,Rc,Rd,k,l,Delta,S00cd,Sklcd)
            call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
            
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)
            call G00ij(k,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_00kl)  
            
            call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)  
            call Gi0j0(i,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i00l)  
            call Gi0j0(i,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i0k0)
            call G0ijk(i,k,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_i0kl)                      
           
            S = Sklcd*( Si0ab*F_0    + S00ab*G_i000 ) + &
                Sk0cd*( Si0ab*G_000l + S00ab*G_i00l ) + &
                S0lcd*( Si0ab*G_00k0 + S00ab*G_i0k0 ) + &
                S00cd*( Si0ab*G_00kl + S00ab*G_i0kl ) 
                
            BNabcd = 2.06409820372476731d0*(a*d*c)**1.25*b**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
        
        if( sa == 0 ) then 
            !------<sp|1/rC|pp>------!
            do j = 1,3
               if( Lb(j) == 1 ) exit
            end do
            
            do k = 1,3
               if( Lc(k) == 1 ) exit
            end do
            
            do l = 1,3
               if( Ld(l) == 1 ) exit
            end do

            call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2) ; call F3(w,F_3)
            call Sab00(a,b,K1,S00ab)
            call Sab00(c,d,K2,S00cd)  
            call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
            call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
            call Sabij(c,d,Rc,Rd,k,l,Delta,S00cd,Sklcd)
            call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
            
            call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
            call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)
            call G00ij(k,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_00kl)  
            
            call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00) 
            call Gi0j0(j,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0j0l)  
            call Gi0j0(j,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0jk0)
            call G0ijk(j,k,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_0jkl)                      
           
            S = Sklcd*( S0jab*F_0    + S00ab*G_0j00 ) + &
                Sk0cd*( S0jab*G_000l + S00ab*G_0j0l ) + &
                S0lcd*( S0jab*G_00k0 + S00ab*G_0jk0 ) + &
                S00cd*( S0jab*G_00kl + S00ab*G_0jkl ) 
                
            BNabcd = 2.06409820372476731d0*(d*b*c)**1.25*a**0.75
            S = S*BNabcd*LAMBDA
            
            return
        end if
       
    end if
    !------------------------- <pp|1/rC|pp> ----------------------------
    if( ijkl == 4 ) then
          
          call F0(w,F_0) ; call F1(w,F_1) ; call F2(w,F_2)
          call F3(w,F_3) ; call F4(w,F_4)
          
          do i = 1,3
             if( La(i) == 1 ) exit
          end do
            
          do j = 1,3
             if( Lb(j) == 1 ) exit
          end do 
            
          do k = 1,3
             if( Lc(k) == 1 ) exit
          end do
            
          do l = 1,3
             if( Ld(l) == 1 ) exit
          end do   
              
          call Sab00(a,b,K1,S00ab)
          call Sab00(c,d,K2,S00cd)  
          call Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)
          call Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)
          call Sabij(a,b,Ra,Rb,i,j,Delta,S00ab,Sijab)          
          call Sabi0(c,d,Rc,Rd,k,S00cd,Sk0cd)
          call Sab0j(c,d,Rc,Rd,l,S00cd,S0lcd)
          call Sabij(c,d,Rc,Rd,k,l,Delta,S00cd,Sklcd)  
          
          call Gi000(j,Rp,Rq,S2,S4,F_1,G_0j00)
          call Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)
          call Gij00(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_ij00)
          
          call G00i0(l,Rp,Rq,S1,S4,F_1,G_000l)
          call Gi0j0(j,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0j0l)
          call Gi0j0(i,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i00l)
          call Gijk0(i,j,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_ij0l)

          call G00i0(k,Rp,Rq,S1,S4,F_1,G_00k0)  
          call Gi0j0(j,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_0jk0)  
          call Gi0j0(i,k,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i0k0)
          call Gijk0(i,j,k,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_ijk0)
          
          call G00ij(k,l,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_00kl)   
          call G0ijk(j,k,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_0jkl)
          call G0ijk(i,k,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_i0kl)
          call Gijkl(i,j,k,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,F_4,G_ijkl)
                  
          S = Sijab*(Sklcd*F_0    + Sk0cd*G_000l + S0lcd*G_00k0 + S00cd*G_00kl) + &
              Si0ab*(Sklcd*G_0j00 + Sk0cd*G_0j0l + S0lcd*G_0jk0 + S00cd*G_0jkl) + &
              S0jab*(Sklcd*G_i000 + Sk0cd*G_i00l + S0lcd*G_i0k0 + S00cd*G_i0kl) + &
              S00ab*(Sklcd*G_ij00 + Sk0cd*G_ij0l + S0lcd*G_ijk0 + S00cd*G_ijkl)
                      
            BNabcd = 4.12819640744953462d0*(a*b*c*d)**1.25
            S = S*BNabcd*LAMBDA
          
          return
    end if
    !-------------------------------------------------------------------
    
end subroutine E_E_Repulsion_sp
!=======================================================================
subroutine Sab00(a,b,K,S00ab)

    implicit none
    real(kind=8),parameter::Pi = 3.141592653589793d0
    real(kind=8) a,b,K,S00ab
    
    S00ab = (Pi/(a+b))**(1.5d0)*K
    
end subroutine Sab00
!=======================================================================
subroutine Sabi0(a,b,Ra,Rb,i,S00ab,Si0ab)

    implicit none
    integer i
    real(kind=8) a,b,S00ab,Si0ab
    real(kind=8) Ra(3),Rb(3)
    
    Si0ab = -b/(a+b)*(Ra(i)-Rb(i))*S00ab
    
end subroutine Sabi0
!=======================================================================
subroutine Sab0j(a,b,Ra,Rb,j,S00ab,S0jab)

    implicit none
    integer j
    real(kind=8) a,b,S00ab,S0jab
    real(kind=8) Ra(3),Rb(3)
    
    S0jab = -a/(a+b)*(Rb(j)-Ra(j))*S00ab
    
end subroutine Sab0j
!=======================================================================
subroutine Sabij(a,b,Ra,Rb,i,j,Delta,S00ab,Sijab)

    implicit none
    integer i,j,Delta(3,3)
    real(kind=8) a,b,S00ab,Sijab
    real(kind=8) Ra(3),Rb(3)
    
    Sijab = (Delta(i,j)/(2*(a+b)) + a*b/(a+b)**2*(Ra(i)-Rb(i))*(Rb(j)-Ra(j)))*S00ab
    
end subroutine Sabij
!=======================================================================
subroutine Gi000(i,Rp,Rq,S2,S4,F_1,G_i000)

    implicit none
    integer i
    real(kind=8) Rp(3),Rq(3),S2,S4,F_1,G_i000
    
    G_i000 = -S2/S4*(Rp(i)-Rq(i))*F_1
    
end subroutine Gi000
!=======================================================================
subroutine G00i0(i,Rp,Rq,S1,S4,F_1,G_00i0)

    implicit none
    integer i
    real(kind=8) Rp(3),Rq(3),S1,S4,F_1,G_00i0
    
    G_00i0 = -S1/S4*(Rq(i)-Rp(i))*F_1
    
end subroutine G00i0
!=======================================================================
subroutine Gij00(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_ij00)

    implicit none
    integer i,j,Delta(3,3)
    real(kind=8) Rp(3),Rq(3),S1,S2,S4,F_1,F_2,G_ij00
    
    G_ij00 = S2/S4*(S2/S4*(Rp(i)-Rq(i))*(Rp(j)-Rq(j))*F_2 - &
             Delta(i,j)*0.5d0/S1*F_1)
    
end subroutine Gij00
!=======================================================================
subroutine G00ij(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_00ij)

    implicit none
    integer i,j,Delta(3,3)
    real(kind=8) Rp(3),Rq(3),S1,S2,S4,F_1,F_2,G_00ij
    
    G_00ij = S1/S4*(S1/S4*(Rq(i)-Rp(i))*(Rq(j)-Rp(j))*F_2 - &
             Delta(i,j)*0.5d0/S2*F_1)
    
end subroutine G00ij
!=======================================================================
subroutine Gi0j0(i,j,Delta,Rp,Rq,S1,S2,S4,F_1,F_2,G_i0j0)

    implicit none
    integer i,j,Delta(3,3)
    real(kind=8) Rp(3),Rq(3),S1,S2,S4,F_1,F_2,G_i0j0
    
    G_i0j0 = (S1*S2/S4*(Rp(i)-Rq(i))*(Rq(j)-Rp(j))*F_2 + &
             Delta(i,j)*0.5d0*F_1)/S4
    
end subroutine Gi0j0
!=======================================================================
subroutine Gijk0(i,j,k,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_ijk0)

    implicit none
    integer i,j,k,Delta(3,3)
    real(kind=8) Rp(3),Rq(3),S1,S2,S4,F_2,F_3,G_ijk0
    
    G_ijk0 = (-S2/S4**2)*( S1*S2/S4*(Rp(i)-Rq(i))*(Rp(j)-Rq(j))*(Rq(k)-Rp(k))*F_3 + &
             0.5d0*F_2*( Delta(i,j)*(Rp(k)-Rq(k)) + Delta(i,k)*(Rp(j)-Rq(j)) + &
             Delta(j,k)*(Rp(i)-Rq(i)) ) ) 
    
end subroutine Gijk0
!=======================================================================
subroutine G0ijk(i,j,k,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,G_0ijk)

    implicit none
    integer i,j,k,Delta(3,3)
    real(kind=8) Rp(3),Rq(3),S1,S2,S4,F_2,F_3,G_0ijk
    
    G_0ijk = (-S1/S4**2)*( S1*S2/S4*(Rp(i)-Rq(i))*(Rq(j)-Rp(j))*(Rq(k)-Rp(k))*F_3 - &
             0.5d0*F_2*( Delta(i,j)*(Rp(k)-Rq(k)) + Delta(i,k)*(Rp(j)-Rq(j)) + &
             Delta(j,k)*(Rp(i)-Rq(i)) ) ) 
    
end subroutine G0ijk
!=======================================================================
subroutine Gijkl(i,j,k,l,Delta,Rp,Rq,S1,S2,S4,F_2,F_3,F_4,G_ijkl)

    implicit none
    integer i,j,k,l,Delta(3,3)
    real(kind=8) Rp(3),Rq(3),S1,S2,S4,F_2,F_3,F_4,G_ijkl
    real(kind=8) pqi,pqj,pqk,pql
    
    pqi = Rp(i)-Rq(i)
    pqj = Rp(j)-Rq(j)
    pqk = Rp(k)-Rq(k)
    pql = Rp(l)-Rq(l)
    
    G_ijkl = ( (S1*S2/S4)**2*(pqi)*(pqj)*(pqk)*(pql)*F_4 - &
             0.5d0*F_3*S1*S2/S4*( Delta(i,j)*(pqk)*(pql) +  &
             Delta(i,k)*(pqj)*(pql) + Delta(i,l)*(pqj)*(pqk) + &
             Delta(j,k)*(pqi)*(pql) + Delta(j,l)*(pqi)*(pqk) + &
             Delta(k,l)*(pqi)*(pqj) ) + 0.25d0*F_2*( Delta(i,j)*Delta(k,l) + &
             Delta(i,k)*Delta(j,l) + Delta(i,l)*Delta(j,k) ) )/S4**2 
   
end subroutine Gijkl
!=======================================================================
