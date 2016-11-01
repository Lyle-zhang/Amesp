!=======================================================================
!                          Single Ele Int                              !
!=======================================================================
subroutine Overlap_Kinetic_Int(Overlap,Kinetic,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)

    implicit none
    integer N,i,j,k,l,bas_num(N),NAtom(N)
    integer l_Basis(3,N),La(3),Lb(3)
    real(kind=8) Overlap(N,N),Kinetic(N,N),Rt(3,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) STemp0,STemp1,STemp2,BNa,BNb
    real(kind=8) KTemp0,KTemp1,KTemp2,a,b
    
    do i = 1,N
       do j = i,N
          STemp1 = 0.0d0
          KTemp1 = 0.0d0
          do k = 1,bas_num(i)
             STemp2 = 0.0d0
             KTemp2 = 0.0d0
             do l = 1,bas_num(j)
                a = a_Basis(k,i,NAtom(i))
                b = a_Basis(l,j,NAtom(j))
                La(:) = l_Basis(:,i)
                Lb(:) = l_Basis(:,j)
                call Bas_Normal(a,La(1),La(2),La(3),BNa)
                call Bas_Normal(b,Lb(1),Lb(2),Lb(3),BNb)
                call Overlap_Kinetic_GE(a_Basis(k,i,NAtom(i)),a_Basis(l,j,NAtom(j)),&
                        l_Basis(:,i),l_Basis(:,j),Rt(:,i),Rt(:,j),BNa,BNb,STemp0,KTemp0,"Y")
                STemp2 = STemp2 + STemp0*c_Basis(k,i,NAtom(i))*c_Basis(l,j,NAtom(j))
                KTemp2 = KTemp2 + KTemp0*c_Basis(k,i,NAtom(i))*c_Basis(l,j,NAtom(j))
             end do
             STemp1 = STemp1 + STemp2
             KTemp1 = KTemp1 + KTemp2
          end do
          Overlap(j,i) = STemp1   !Overlap Matrix
          Overlap(i,j) = STemp1
          Kinetic(j,i) = KTemp1   !Kinetic Matrix
          Kinetic(i,j) = KTemp1
       end do
    end do

end subroutine Overlap_Kinetic_Int
!=======================================================================
subroutine N_E_Attractive_Int(N_E_Attract,M,N,a_Basis,c_Basis,l_Basis,bas_num,Rn,Rt,MAtom,NAtom)

    implicit none
    integer M,N,i,j,k,l,bas_num(N),MAtom(M),NAtom(N),l_Basis(3,N)
    real(kind=8) N_E_Attract(N,N),Rn(3,M),Rt(3,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) Temp0,Temp1,Temp2
        
    do i = 1,N
      do j = i,N
        Temp1 = 0.0d0
        do k = 1,bas_num(i)
           Temp2 = 0.0d0
           do l = 1,bas_num(j)
                !call N_E_Attractive_GE(a,b,La,Lb,Ra,Rb,Rn,M,MAtom,S)
                call N_E_Attractive_GE(a_Basis(k,i,NAtom(i)),a_Basis(l,j,NAtom(j)), &
                          l_Basis(:,i),l_Basis(:,j),Rt(:,i),Rt(:,j),Rn,M,MAtom,Temp0)                   
                Temp2 = Temp2 + Temp0*c_Basis(k,i,NAtom(i))*c_Basis(l,j,NAtom(j))
           end do
           Temp1 = Temp1 + Temp2
        end do 
        N_E_Attract(j,i) = -Temp1      !N_E_Attract Matrix 
        N_E_Attract(i,j) = -Temp1
      end do 
    end do 

end subroutine N_E_Attractive_Int
!=======================================================================
subroutine Dipole_Int(Dipole,Overlap,N,a_Basis,c_Basis,l_Basis,bas_num,Rt,NAtom)

    implicit none
    integer N,i,j,k,l,bas_num(N),NAtom(N)
    integer l_Basis(3,N),La(3),Lb(3)
    real(kind=8) Overlap(N,N),Dipole(N,N,3),Rt(3,N)
    real(kind=8) c_Basis(6,N,36),a_Basis(6,N,36)
    real(kind=8) STemp0,STemp1,STemp2
    real(kind=8) KTemp0,BNa,BNb,a,b
    
    do i = 1,N
       do j = i,N
          STemp1 = 0.0d0
          do k = 1,bas_num(i)
             STemp2 = 0.0d0
             do l = 1,bas_num(j)
                a = a_Basis(k,i,NAtom(i))
                b = a_Basis(l,j,NAtom(j))
                La(:) = l_Basis(:,i)
                Lb(:) = l_Basis(:,j)
                call Bas_Normal(a,La(1),La(2),La(3),BNa)
                call Bas_Normal(b,Lb(1),Lb(2),Lb(3),BNb)
                Lb(1) = Lb(1) + 1
                call Overlap_Kinetic_GE(a_Basis(k,i,NAtom(i)),a_Basis(l,j,NAtom(j)),&
                                       La,Lb,Rt(:,i),Rt(:,j),BNa,BNb,STemp0,KTemp0,"N")
                STemp2 = STemp2 + STemp0*c_Basis(k,i,NAtom(i))*c_Basis(l,j,NAtom(j))
             end do
             STemp1 = STemp1 + STemp2
          end do
          Dipole(j,i,1) = STemp1 + Rt(1,j)*Overlap(j,i)   ! Dipole(x)
          Dipole(i,j,1) = Dipole(j,i,1)
       end do
    end do
    !-------------------------------------------------------------------
    do i = 1,N
       do j = i,N
          STemp1 = 0.0d0
          do k = 1,bas_num(i)
             STemp2 = 0.0d0
             do l = 1,bas_num(j)
                a = a_Basis(k,i,NAtom(i))
                b = a_Basis(l,j,NAtom(j))
                La(:) = l_Basis(:,i)
                Lb(:) = l_Basis(:,j)
                call Bas_Normal(a,La(1),La(2),La(3),BNa)
                call Bas_Normal(b,Lb(1),Lb(2),Lb(3),BNb)
                Lb(2) = Lb(2) + 1
                call Overlap_Kinetic_GE(a_Basis(k,i,NAtom(i)),a_Basis(l,j,NAtom(j)),&
                                       La,Lb,Rt(:,i),Rt(:,j),BNa,BNb,STemp0,KTemp0,"N")
                STemp2 = STemp2 + STemp0*c_Basis(k,i,NAtom(i))*c_Basis(l,j,NAtom(j))
             end do
             STemp1 = STemp1 + STemp2
          end do
          Dipole(j,i,2) = STemp1 + Rt(2,j)*Overlap(j,i)    ! Dipole(y)
          Dipole(i,j,2) = Dipole(j,i,2)
       end do
    end do
    !-------------------------------------------------------------------
    do i = 1,N
       do j = i,N
          STemp1 = 0.0d0
          do k = 1,bas_num(i)
             STemp2 = 0.0d0
             do l = 1,bas_num(j)
                a = a_Basis(k,i,NAtom(i))
                b = a_Basis(l,j,NAtom(j))
                La(:) = l_Basis(:,i)
                Lb(:) = l_Basis(:,j)
                call Bas_Normal(a,La(1),La(2),La(3),BNa)
                call Bas_Normal(b,Lb(1),Lb(2),Lb(3),BNb)
                Lb(3) = Lb(3) + 1
                call Overlap_Kinetic_GE(a_Basis(k,i,NAtom(i)),a_Basis(l,j,NAtom(j)),&
                                       La,Lb,Rt(:,i),Rt(:,j),BNa,BNb,STemp0,KTemp0,"N")
                STemp2 = STemp2 + STemp0*c_Basis(k,i,NAtom(i))*c_Basis(l,j,NAtom(j))
             end do
             STemp1 = STemp1 + STemp2
          end do
          Dipole(j,i,3) = STemp1 + Rt(3,j)*Overlap(j,i)   ! Dipole(z)
          Dipole(i,j,3) = Dipole(j,i,3)
       end do
    end do
    !-------------------------------------------------------------------

end subroutine Dipole_Int
!=======================================================================
!                                                                      !
!=======================================================================
subroutine Overlap_Kinetic_GE(a,b,La,Lb,Ra,Rb,BNa,BNb,S,Ki,YN)

    implicit none
    integer i,i1,j1
    integer La(3),Lb(3),Ca,Cb,IT
    real(kind=8),parameter::Pi = 3.141592653589793d0
    real(kind=8) a,b,K,AB,S,Ki,Ix,Iy,Iz,KIx,KIy,KIz
    real(kind=8) I2x,I2y,I2z,I_2x,I_2y,I_2z
    real(kind=8) BNa,BNb,fx(0:La(1)+Lb(1)+2),fy(0:La(2)+Lb(2)+2),fz(0:La(3)+Lb(3)+2)
    real(kind=8) Ra(3),Rb(3),Rp(3),PA(3),PB(3),Temp
    character(len=1) YN
    
    call GTO_GE(a,b,Ra,Rb,K,AB,Rp,PA,PB)

    !------------------------------<0|0>--------------------------------
    do i = 0,La(1)+Lb(1)
       if( mod(i,2) == 1 ) cycle
       fx(i) = 0.0d0
       do i1 = 0,La(1)
          do j1 = 0,Lb(1)
             if ( i1+j1 == i ) then
                call Combination(La(1),i1,Ca)
                call Combination(Lb(1),j1,Cb)
                fx(i) = fx(i) + Ca*Cb*PA(1)**(La(1)-i1)*PB(1)**(Lb(1)-j1)
             end if
          end do
       end do
    end do
    Ix = 0.0d0
    do i = 0,(La(1)+Lb(1))/2
        call Double_Factorial(2*i-1,IT)
        Ix = Ix + fx(i*2)*IT/(2.0d0*(a+b))**i
    end do

    
    do i = 0,La(2)+Lb(2)
       if( mod(i,2) == 1 ) cycle
       fy(i) = 0.0d0
       do i1 = 0,La(2)
          do j1 = 0,Lb(2)
             if ( i1+j1 == i ) then
                call Combination(La(2),i1,Ca)
                call Combination(Lb(2),j1,Cb)
                fy(i) = fy(i) + Ca*Cb*PA(2)**(La(2)-i1)*PB(2)**(Lb(2)-j1)
             end if
          end do
       end do
    end do    
    Iy = 0.0d0
    do i = 0,(La(2)+Lb(2))/2
        call Double_Factorial(2*i-1,IT)
        Iy = Iy + fy(i*2)*IT/(2.0d0*(a+b))**i
    end do

        
    do i = 0,La(3)+Lb(3)
       if( mod(i,2) == 1 ) cycle
       fz(i) = 0.0d0
       do i1 = 0,La(3)
          do j1 = 0,Lb(3)
             if ( i1+j1 == i ) then
                call Combination(La(3),i1,Ca)
                call Combination(Lb(3),j1,Cb)
                fz(i) = fz(i) + Ca*Cb*PA(3)**(La(3)-i1)*PB(3)**(Lb(3)-j1)
             end if
          end do
       end do
    end do
    Iz = 0.0d0
    do i = 0,(La(3)+Lb(3))/2
        call Double_Factorial(2*i-1,IT)
        Iz = Iz + fz(i*2)*IT/(2.0d0*(a+b))**i
    end do
    
    Temp = (Pi/(a+b))**(1.5d0)*BNa*BNb*K
    S = Temp*Ix*Iy*Iz 
    if( YN /= "Y" ) return
    !------------------------------<0|+2>-------------------------------
    do i = 0,La(1)+Lb(1)+2
       if( mod(i,2) == 1 ) cycle
       fx(i) = 0.0d0
       do i1 = 0,La(1)
          do j1 = 0,Lb(1)+2
             if ( i1+j1 == i ) then
                call Combination(La(1),i1,Ca)
                call Combination(Lb(1)+2,j1,Cb)
                fx(i) = fx(i) + Ca*Cb*PA(1)**(La(1)-i1)*PB(1)**(Lb(1)+2-j1)
             end if
          end do
       end do
    end do
    I2x = 0.0d0
    do i = 0,(La(1)+Lb(1))/2+1
        call Double_Factorial(2*i-1,IT)
        I2x = I2x + fx(i*2)*IT/(2.0d0*(a+b))**i
    end do

    
    do i = 0,La(2)+Lb(2)+2
       if( mod(i,2) == 1 ) cycle
       fy(i) = 0.0d0
       do i1 = 0,La(2)
          do j1 = 0,Lb(2)+2
             if ( i1+j1 == i ) then
                call Combination(La(2),i1,Ca)
                call Combination(Lb(2)+2,j1,Cb)
                fy(i) = fy(i) + Ca*Cb*PA(2)**(La(2)-i1)*PB(2)**(Lb(2)+2-j1)
             end if
          end do
       end do
    end do    
    I2y = 0.0d0
    do i = 0,(La(2)+Lb(2))/2+1
        call Double_Factorial(2*i-1,IT)
        I2y = I2y + fy(i*2)*IT/(2.0d0*(a+b))**i
    end do

        
    do i = 0,La(3)+Lb(3)+2
       if( mod(i,2) == 1 ) cycle
       fz(i) = 0.0d0
       do i1 = 0,La(3)
          do j1 = 0,Lb(3)+2
             if ( i1+j1 == i ) then
                call Combination(La(3),i1,Ca)
                call Combination(Lb(3)+2,j1,Cb)
                fz(i) = fz(i) + Ca*Cb*PA(3)**(La(3)-i1)*PB(3)**(Lb(3)+2-j1)
             end if
          end do
       end do
    end do
    I2z = 0.0d0
    do i = 0,(La(3)+Lb(3))/2+1
        call Double_Factorial(2*i-1,IT)
        I2z = I2z + fz(i*2)*IT/(2.0d0*(a+b))**i
    end do
    !------------------------------<0|-2>-------------------------------
    if ( Lb(1) < 2 ) then
       I_2x = 0.0d0
    else   
    do i = 0,La(1)+Lb(1)-2
       if( mod(i,2) == 1 ) cycle
       fx(i) = 0.0d0
       do i1 = 0,La(1)
          do j1 = 0,Lb(1)-2
             if ( i1+j1 == i ) then
                call Combination(La(1),i1,Ca)
                call Combination(Lb(1)-2,j1,Cb)
                fx(i) = fx(i) + Ca*Cb*PA(1)**(La(1)-i1)*PB(1)**(Lb(1)-2-j1)
             end if
          end do
       end do
    end do
    I_2x = 0.0d0
    do i = 0,(La(1)+Lb(1))/2-1
        call Double_Factorial(2*i-1,IT)
        I_2x = I_2x + fx(i*2)*IT/(2.0d0*(a+b))**i
    end do
    end if

    if ( Lb(2) < 2 ) then
       I_2y = 0.0d0
    else     
    do i = 0,La(2)+Lb(2)-2
       if( mod(i,2) == 1 ) cycle
       fy(i) = 0.0d0
       do i1 = 0,La(2)
          do j1 = 0,Lb(2)-2
             if ( i1+j1 == i ) then
                call Combination(La(2),i1,Ca)
                call Combination(Lb(2)-2,j1,Cb)
                fy(i) = fy(i) + Ca*Cb*PA(2)**(La(2)-i1)*PB(2)**(Lb(2)-2-j1)
             end if
          end do
       end do
    end do    
    I_2y = 0.0d0
    do i = 0,(La(2)+Lb(2))/2-1
        call Double_Factorial(2*i-1,IT)
        I_2y = I_2y + fy(i*2)*IT/(2.0d0*(a+b))**i
    end do
    end if

    if ( Lb(3) < 2 ) then
       I_2z = 0.0d0
    else         
    do i = 0,La(3)+Lb(3)-2
       if( mod(i,2) == 1 ) cycle
       fz(i) = 0.0d0
       do i1 = 0,La(3)
          do j1 = 0,Lb(3)-2
             if ( i1+j1 == i ) then
                call Combination(La(3),i1,Ca)
                call Combination(Lb(3)-2,j1,Cb)
                fz(i) = fz(i) + Ca*Cb*PA(3)**(La(3)-i1)*PB(3)**(Lb(3)-2-j1)
             end if
          end do
       end do
    end do
    I_2z = 0.0d0
    do i = 0,(La(3)+Lb(3))/2-1
        call Double_Factorial(2*i-1,IT)
        I_2z = I_2z + fz(i*2)*IT/(2.0d0*(a+b))**i
    end do
    end if
    !-------------------------------------------------------------------
    I2x = Temp*I2x*Iy*Iz
    I2y = Temp*Ix*I2y*Iz
    I2z = Temp*Ix*Iy*I2z
    
    I_2x = Temp*I_2x*Iy*Iz
    I_2y = Temp*Ix*I_2y*Iz
    I_2z = Temp*Ix*Iy*I_2z

    KIx = -0.5d0*Lb(1)*(Lb(1)-1)*I_2x + b*(2*Lb(1)+1)*S - 2.0d0*b**2*I2x
    KIy = -0.5d0*Lb(2)*(Lb(2)-1)*I_2y + b*(2*Lb(2)+1)*S - 2.0d0*b**2*I2y
    KIz = -0.5d0*Lb(3)*(Lb(3)-1)*I_2z + b*(2*Lb(3)+1)*S - 2.0d0*b**2*I2z
    !-------------------------------------------------------------------
   
    Ki = KIx + KIy + KIz

end subroutine Overlap_Kinetic_GE
!=======================================================================
subroutine N_E_Attractive_GE(a,b,La,Lb,Ra,Rb,Rn,M,MAtom,S)

    implicit none
    integer ,parameter :: LN = 4 
    integer La(3),Lb(3),ii,i,j,k,M,MAtom(M)
    integer L_max,Lx,Ly,Lz
    real(kind=8),parameter::Pi = 3.141592653589793d0
    real(kind=8) a,b,AB,BNa,BNb,S
    real(kind=8) K0,Rp(3),PA(3),PB(3),PC(3),T
    real(kind=8) RPC,w,F(0:LN),Gx(0:LN),Gy(0:LN),Gz(0:LN)
    real(kind=8) Rn(3,M),Ra(3),Rb(3),Rc(3)
    real(kind=8) Temp1,Temp2,Temp3
    
    call GTO_GE(a,b,Ra,Rb,K0,AB,Rp,PA,PB)
    
    call Bas_Normal(a,La(1),La(2),La(3),BNa)
    call Bas_Normal(b,Lb(1),Lb(2),Lb(3),BNb)
    
    L_max = sum(La(:)) + sum(Lb(:))
    Lx = La(1) + Lb(1) ; Ly = La(2) + Lb(2) ; Lz = La(3) + Lb(3) 
    T = 1.0d0/(a+b)
    
    S = 0.0d0
    do ii = 1,M
    
        Rc = Rn(:,ii)
        PC = - Rc + Rp
        RPC = (Rp(1)-Rc(1))**2+(Rp(2)-Rc(2))**2+(Rp(3)-Rc(3))**2
        w = RPC*(a+b)
        
        do i = 0,L_max
            call Boys(i,w,F(i))
        end do
        
        call G_Func(la(1),lb(1),PA(1),PB(1),PC(1),T,Gx(0:Lx))
        call G_Func(la(2),lb(2),PA(2),PB(2),PC(2),T,Gy(0:Ly))      
        call G_Func(la(3),lb(3),PA(3),PB(3),PC(3),T,Gz(0:Lz))
        
        !---------------------------------------------------------------
        Temp1 = 0.0d0
        do i = 0,Lx
            Temp2 = 0.0d0
            do j = 0,Ly
                Temp3 = 0.0d0
                do k = 0,Lz
                    Temp3 = Temp3 + Gx(i)*Gy(j)*Gz(k)*F(i+j+k)
                end do
                Temp2 = Temp2 + Temp3
            end do
            Temp1 = Temp1 + Temp2
        end do
        !---------------------------------------------------------------
        S = S + Temp1*MAtom(ii)
        
    end do
    
    S = S*K0*2.0d0*Pi/(a+b)*BNa*BNb

end subroutine N_E_Attractive_GE
!=======================================================================
subroutine G_Func(l1,l2,a1,b1,c,T,G)

    implicit none
    integer l1,l2,L
    real(kind=8) a,b,a1,b1,c,T,G(0:l1+l2)
    L = l1 + l2
    
    if( l1 < l2 ) then
        a = b1
        b = a1
    else
        a = a1
        b = b1
    end if
    
    select case(L)
    
        case(0)
            G(0) = 1.0d0
        case(1)
            G(0) = a
            G(1) = -c
        case(2)
            if(l1 == l2) then
                G(0) = a*b + T*0.5
                G(1) = -(a+b)*c - T*0.5d0
                G(2) = c*c
            else
                G(0) = a*a + T*0.5
                G(1) = -2*a*c - T*0.5d0
                G(2) = c*c
            end if
        case(3)
            G(0) = a*a*b + (2*a + b)*T*0.5
            G(1) = -(a*a + 2*a*b)*c - (2*a + b)*T*0.5d0 - 1.5d0*T*c
            G(2) = (2*a+b)*c*c + 1.5d0*T*c
            G(3) = -c*c*c
        case(4)
            G(0) = a*a*b*b + (a*a + 4*a*b + b*b)*T*0.5d0 + 0.75d0*T*T
            G(1) = -2*(a*a*b + a*b*b)*c - (a*a + 4*a*b + b*b)*T*0.5d0 - 3*(a + b)*c*T - 1.5d0*T*T
            G(2) = (a*a + 4*a*b + b*b)*c*c + 3*c*c*T + 3*(a+b)*c*T + 0.75d0*T*T
            G(3) = -2*(a+b)*c*c*c - 3*c*c*T
            G(4) = c**4
    end select
    
    return
    
end subroutine G_Func
!=======================================================================

