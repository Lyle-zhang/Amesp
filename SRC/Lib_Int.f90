!=======================================================================
!                               Int Dep                                !
!=======================================================================
subroutine GTO_s(a,b,Ra,Rb,K,AB,Rp) 

    implicit none
    integer i
    real(kind=8) a,b,K
    real(kind=8) Ra(3),Rb(3),Rp(3)
    real(kind=8) AB
    
    AB = (Ra(1)-Rb(1))**2+(Ra(2)-Rb(2))**2+(Ra(3)-Rb(3))**2
    K = dexp(-a*b/(a+b)*AB)
    
    do i = 1,3
         Rp(i) = (Rb(i)-Ra(i))*(b/(a+b)) + Ra(i)
    end do
    
end subroutine GTO_s
!=======================================================================
subroutine GTO_GE(a,b,Ra,Rb,K,AB,Rp,PA,PB) 

    implicit none
    integer i
    real(kind=8) a,b,K
    real(kind=8) Ra(3),Rb(3),Rp(3),PA(3),PB(3)
    real(kind=8) AB,cosi(3)
    
    AB = (Ra(1)-Rb(1))**2+(Ra(2)-Rb(2))**2+(Ra(3)-Rb(3))**2
    
    if( AB < 1.0d-12 ) then
       K = 1.0d0
       PA = 0.0d0
       PB = 0.0d0
       Rp = Ra
       return
    end if
    
    K = dexp(-a*b/(a+b)*AB)
    
    do i = 1,3
         Rp(i) = (Ra(i)*a + Rb(i)*b)/(a+b) 
         cosi(i) = (Rb(i)-Ra(i))
         PA(i) = (b/(a+b))*cosi(i)
         PB(i) = -(a/(a+b))*cosi(i)
    end do   
    
end subroutine GTO_GE
!=======================================================================
subroutine Bas_Normal(a,l,m,n,BN)

    implicit none
    integer l,m,n,lmn
    integer L2,M2,N2
    real(kind=8),parameter::Pi = 3.141592653589793d0
    real(kind=8) a,BN
    
    lmn = l + m + n
    
    if( lmn == 0 ) then
        BN = 0.7127054703549902d0*a**0.75
        return
    end if
    
    if( lmn == 1) then
        BN = 1.4254109407099803d0*a**1.25
        return
    end if
    
    if( lmn == 2) then
        if( l == 2 .or. m == 2 .or. n == 2 ) then
            BN = 1.64592278064948966d0*a**1.75
            return
        else
            BN = 2.85082188141996064d0*a**1.75
            return
        end if
    end if
    
    if( lmn > 2 ) then
    
        call Double_Factorial(2*l-1,L2)
        call Double_Factorial(2*m-1,M2)
        call Double_Factorial(2*n-1,N2)
    
        BN = (2.0d0*a/Pi)**(0.75)*dsqrt((4.0d0*a)**(l+m+n)/(real(L2*M2*N2,8)))
        return
        
    end if
    
end subroutine Bas_Normal
!=======================================================================
subroutine Double_Factorial(N,A)

    implicit none    
    integer A
    integer i,N,Z
    
    A = 1
    
    if(N > 16) then
        write(*,*) "N must < 16"
        stop
    end if
    if(N == 0) then 
        A = 1
    else 
        Z = mod(N,2)
        select case(Z)
          case(1)
             do i = 1,N,2
                A = A*i
             end do
          case(0)
             do i = 2,N,2
                A = A*i
             end do
        end select
    end if
    
 end subroutine Double_Factorial
!=======================================================================
subroutine Combination(n,m,P)
 
    implicit none
    integer A1,A2,A3,P
    integer m,n
    call Factorial(N,A1)
    call Factorial(N-M,A2)
    call Factorial(M,A3)
    P = A1/(A2*A3)
    
end subroutine
!=======================================================================
subroutine Factorial(N,A)

    implicit none   
    integer i,N,A  
    A = 1
    if(N==0) return
    if(N>14) then
        write(*,*) "N must < 14"
        stop
    end if
    
    do i = 1,N
         A = A*i
    end do

end subroutine Factorial
!=======================================================================
!                               Boys                                   !
!=======================================================================
subroutine Boys(m,w,F)
   
    implicit none
    integer m
    real(kind=8)  w,F
    
    if ( m == 0 ) then
       call F0(w,F)
       return
    end if
    
    if ( m == 1 ) then
       call F1(w,F)
       return
    end if    
    
    if ( m == 2 ) then
       call F2(w,F)
       return
    end if 
    
    if ( m == 3 ) then
       call F3(w,F)
       return
    end if 
    
    if ( m == 4 ) then
       call F4(w,F)
       return
    end if 
          
end subroutine Boys
!=======================================================================
subroutine F0(w,F)

    implicit none
    real(kind=8),parameter:: w0 = 16.3578d0
    
    real(kind=8),parameter:: a1 = 0.213271302431420d0
    real(kind=8),parameter:: a2 = 0.629344460255614d-1
    real(kind=8),parameter:: a3 = 0.769838037756759d-2
    real(kind=8),parameter:: a4 = 0.758433197127160d-3
    real(kind=8),parameter:: a5 = 0.564691197633667d-4
    
    real(kind=8),parameter:: b1 = 0.879937801660182d0
    real(kind=8),parameter:: b2 = 0.338450368470103d0
    real(kind=8),parameter:: b3 = 0.738522953299624d-1
    real(kind=8),parameter:: b4 = 0.101431553402629d-1
    real(kind=8),parameter:: b5 = 0.955528842975585d-3
    real(kind=8),parameter:: b6 = 0.720266520392572d-4
    
    real(kind=8) w,F
    real(kind=8) Tempa,Tempb
    
    if( w > w0 ) then
        F = 0.88622692545275801d0/dsqrt(w)
        return
    else if( w > 1d-8 .and. w <= w0 ) then
        Tempa = 1.0d0+(a1+(a2+(a3+(a4+a5*w)*w)*w)*w)*w
        Tempb = 1.0d0+(b1+(b2+(b3+(b4+(b5+b6*w)*w)*w)*w)*w)*w
        F = dsqrt(Tempa/Tempb)
        return
    else if( w < 1d-8 ) then
        F = 1.0d0
        return
    end if
          
end subroutine F0
!=======================================================================
subroutine F1(w,F)

    implicit none
    real(kind=8),parameter:: w0 = 17.4646d0
    
    real(kind=8),parameter:: a0 = 0.480749856769136d0
    real(kind=8),parameter:: a1 = 0.295195994716045d-1
    real(kind=8),parameter:: a2 = 0.128790985465415d-1
    real(kind=8),parameter:: a3 = 0.998165499553218d-3
    real(kind=8),parameter:: a4 = 0.970927983276419d-4
    real(kind=8),parameter:: a5 = 0.493839847029699d-5
    
    real(kind=8),parameter:: b1 = 0.461403194579124d0
    real(kind=8),parameter:: b2 = 0.108494164372449d0
    real(kind=8),parameter:: b3 = 0.171462934845042d-1
    real(kind=8),parameter:: b4 = 0.196918657845508d-2
    real(kind=8),parameter:: b5 = 0.160138863265245d-3
    real(kind=8),parameter:: b6 = 0.857708713007233d-5
    
    real(kind=8) w,F
    real(kind=8) Tempa,Tempb
    
    if( w > w0 ) then
        F = 0.443113462726379007d0/(w**(1.5d0))
        return
    else if( w > 1d-8 .and. w <= w0 ) then
        Tempa = a0+(a1+(a2+(a3+(a4+a5*w)*w)*w)*w)*w
        Tempb = 1.0d0+(b1+(b2+(b3+(b4+(b5+b6*w)*w)*w)*w)*w)*w
        F = (Tempa/Tempb)**1.5d0
        return
    else if( w < 1d-8 ) then
        F = 0.33333333333333333d0
        return
    end if
          
end subroutine F1
!=======================================================================
subroutine F2(w,F)

    implicit none
    real(kind=8),parameter:: w0 = 15.2368d0
    
    real(kind=8),parameter:: a0 = 0.5253055608807534d0
    real(kind=8),parameter:: a1 = -0.575763488635418d-2
    real(kind=8),parameter:: a2 = 0.731474973333076d-2
    real(kind=8),parameter:: a3 = 0.251276149443393d-3
    real(kind=8),parameter:: a4 = 0.264336244559094d-4
    
    real(kind=8),parameter:: b1 = 0.274754154712841d0
    real(kind=8),parameter:: b2 = 0.425364830353043d-1
    real(kind=8),parameter:: b3 = 0.493902790955943d-2
    real(kind=8),parameter:: b4 = 0.437251500927601d-3
    real(kind=8),parameter:: b5 = 0.288914662393981d-4
    
    real(kind=8) w,F
    real(kind=8) Tempa,Tempb
    
    if( w > w0 ) then
        F = 0.6646701940895684d0/(w**(2.5d0))
        return
    else if( w > 1d-8 .and. w <= w0 ) then
        Tempa = a0+(a1+(a2+(a3+a4*w)*w)*w)*w
        Tempb = 1.0d0+(b1+(b2+(b3+(b4+b5*w)*w)*w)*w)*w
        F = (Tempa/Tempb)**2.5d0
        return
    else if( w < 1d-8 ) then
        F = 0.2d0
        return
    end if
          
end subroutine F2
!=======================================================================
subroutine F3(w,F)

    implicit none
    real(kind=8),parameter:: w0 = 16.0419d0
    
    real(kind=8),parameter:: a0 = 0.5735131987446478d0
    real(kind=8),parameter:: a1 = -0.290110430424666d-1
    real(kind=8),parameter:: a2 = 0.561884370781462d-2
    real(kind=8),parameter:: a3 = 0.301628267382713d-4
    real(kind=8),parameter:: a4 = 0.110671035361856d-4
    
    real(kind=8),parameter:: b1 = 0.171637608242892d0
    real(kind=8),parameter:: b2 = 0.187571417256877d-1
    real(kind=8),parameter:: b3 = 0.178536829675118d-2
    real(kind=8),parameter:: b4 = 0.137360778130936d-3
    real(kind=8),parameter:: b5 = 0.791915206883054d-5
    
    real(kind=8) w,F
    real(kind=8) Tempa,Tempb
    
    if( w > w0 ) then
        F = 1.661675485223921d0/(w**(3.5d0))
        return
    else if( w > 1d-8 .and. w <= w0 ) then
        Tempa = a0+(a1+(a2+(a3+a4*w)*w)*w)*w
        Tempb = 1.0d0+(b1+(b2+(b3+(b4+b5*w)*w)*w)*w)*w
        F = (Tempa/Tempb)**3.5d0
        return
    else if( w < 1d-8 ) then
        F = 0.14285714285714288d0
        return
    end if
          
end subroutine F3
!=======================================================================
subroutine F4(w,F)

    implicit none
    real(kind=8),parameter:: w0 = 16.8955d0
    
    real(kind=8),parameter:: a0 = 0.613685849032916d0
    real(kind=8),parameter:: a1 = -0.452693111179624d-1
    real(kind=8),parameter:: a2 = 0.490070062899003d-2
    real(kind=8),parameter:: a3 = -0.561789719979307d-4
    real(kind=8),parameter:: a4 = 0.550814626951998d-5
    
    real(kind=8),parameter:: b1 = 0.108051989937231d0
    real(kind=8),parameter:: b2 = 0.855924943430755d-2
    real(kind=8),parameter:: b3 = 0.724968571389473d-3
    real(kind=8),parameter:: b4 = 0.502338223156067d-4
    real(kind=8),parameter:: b5 = 0.249107837399141d-5
    
    real(kind=8) w,F
    real(kind=8) Tempa,Tempb
    
    if( w > w0 ) then
        F = 5.815864198283723d0/(w**(4.5d0))
        return
    else if( w > 1d-8 .and. w <= w0 ) then
        Tempa = a0+(a1+(a2+(a3+a4*w)*w)*w)*w
        Tempb = 1.0d0+(b1+(b2+(b3+(b4+b5*w)*w)*w)*w)*w
        F = (Tempa/Tempb)**4.5d0
        return
    else if( w < 1d-8 ) then
        F = 0.1111111111111111d0
        return
    end if
          
end subroutine F4
!=======================================================================
