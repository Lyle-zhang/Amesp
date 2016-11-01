subroutine MP2_Close(N,Nocc,En,g_Mo,EI,mp_2)
    !------------- Close Shell,Multiplicity = 1 --------------!
    implicit none
    integer N,Nocc,k,l,a,b,EI(N,N)
    real(kind=8) mp_2
    real(kind=8) En(N)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) Temp1,Temp2,Temp3,Temp4
               
    Temp1 = 0.0d0
    do k = 1,Nocc
      Temp2 = 0.0d0
      do l = 1,Nocc
        Temp3 = 0.0d0
        do a = Nocc+1,N
          Temp4 = 0.0d0
          do b = Nocc+1,N        
             Temp4 = Temp4 + (g_Mo(EI(b,l),EI(a,k))*(2.0d0*g_Mo(EI(l,b),EI(k,a)) &
                     -g_Mo(EI(l,a),EI(k,b))))/(En(k)+En(l)-En(a)-En(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp_2 = Temp1
    
end subroutine MP2_Close
!=======================================================================
subroutine MP3_Close(N,Nocc,En,g_Mo,EI,mp_3)
    !------------- Close Shell,Multiplicity = 1 --------------!
    implicit none
    integer N,Nocc,k,l,m,n0,a,b,c,d,EI(N,N)
    real(kind=8) mp_3,En(N)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) Temp1,Temp2,Temp3,Temp4,Temp5,Temp6
    real(kind=8) Rklab,Rkmac,Rklcd,Rmnab
    real(kind=8) e3_I,e3_II,e3_III
    
    Temp1 = 0.0d0
    do k = 1,Nocc
     Temp2 = 0.0d0
     do l = 1,Nocc
      Temp3 = 0.0d0
      do m = 1,Nocc
       Temp4 = 0.0d0
       do a = Nocc+1,N
        Temp5 = 0.0d0
        do b = Nocc+1,N
         Temp6 = 0.0d0
         do c = Nocc+1,N
            Rklab = 1.0d0/(En(k)+En(l)-En(a)-En(b))
            Rkmac = 1.0d0/(En(k)+En(m)-En(a)-En(c))
            Temp6 = Temp6 + Rklab*Rkmac*((2.0d0*g_Mo(EI(b,l),EI(a,k))-g_Mo(EI(a,l),EI(b,k)))* &
                    (2.0d0*g_Mo(EI(c,m),EI(b,l))-g_Mo(EI(c,b),EI(m,l)))*(2.0d0*g_Mo(EI(c,m),EI(a,k))- &
                    g_Mo(EI(a,m),EI(c,k)))-3.0d0*(g_Mo(EI(a,l),EI(b,k))*g_Mo(EI(c,b),EI(m,l))* &
                    g_Mo(EI(a,m),EI(c,k))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    e3_I = Temp1
    
    Temp1 = 0.0d0
    do k = 1,Nocc
     Temp2 = 0.0d0
     do l = 1,Nocc
      Temp3 = 0.0d0
      do a = Nocc+1,N
       Temp4 = 0.0d0
       do b = Nocc+1,N
        Temp5 = 0.0d0
        do c = Nocc+1,N
         Temp6 = 0.0d0
         do d = Nocc+1,N
            Rklab = 1.0d0/(En(k)+En(l)-En(a)-En(b))
            Rklcd = 1.0d0/(En(k)+En(l)-En(c)-En(d))
            Temp6 = Temp6 + Rklab*Rklcd*g_Mo(EI(b,l),EI(a,k))*g_Mo(EI(d,b),EI(c,a))* &
                    (2.0d0*g_Mo(EI(d,l),EI(c,k))-g_Mo(EI(c,l),EI(d,k)))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    e3_II = Temp1
    
    Temp1 = 0.0d0
    do k = 1,Nocc
     Temp2 = 0.0d0
     do l = 1,Nocc
      Temp3 = 0.0d0
      do m = 1,Nocc
       Temp4 = 0.0d0
       do n0 = 1,Nocc
        Temp5 = 0.0d0
        do a = Nocc+1,N
         Temp6 = 0.0d0
         do b = Nocc+1,N
            Rklab = 1.0d0/(En(k)+En(l)-En(a)-En(b))
            Rmnab = 1.0d0/(En(m)+En(n0)-En(a)-En(b))
            Temp6 = Temp6 + Rklab*Rmnab*g_Mo(EI(b,l),EI(a,k))*g_Mo(EI(n0,l),EI(m,k))* &
                    (2.0d0*g_Mo(EI(b,n0),EI(a,m))-g_Mo(EI(b,m),EI(a,n0)))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    e3_III = Temp1
        
    mp_3 = e3_I + e3_II + e3_III
    
end subroutine MP3_Close
!=======================================================================
!=======================================================================
subroutine MP2_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_2)
    !------------- Open Shell,Multiplicity > 1 --------------!
    implicit none
    integer N,NoccA,NoccB,k,l,a,b,EI(N,N)
    real(kind=8) mp_2,mp2_Temp(6)
    real(kind=8) EnA(N),EnB(N),g_MoAB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) Temp1,Temp2,Temp3,Temp4

    Temp1 = 0.0d0
    do k = 1,NoccB
      Temp2 = 0.0d0
      do l = 1,NoccA
        Temp3 = 0.0d0
        do a = NoccA+1,N
          Temp4 = 0.0d0
          do b = NoccB+1,N        
             Temp4 = Temp4 + g_MoAB(EI(b,k),EI(a,l))*g_MoAB(EI(b,k),EI(a,l))/(EnB(k)+EnA(l)-EnA(a)-EnB(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp2_Temp(1) = Temp1  
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do k = 1,NoccB
      Temp2 = 0.0d0
      do l = 1,NoccB
        Temp3 = 0.0d0
        do a = NoccB+1,N
          Temp4 = 0.0d0
          do b = NoccB+1,N        
             Temp4 = Temp4 + (g_MoB(EI(b,l),EI(a,k))-g_MoB(EI(a,l),EI(b,k)))**2/(EnB(k)+EnB(l)-EnB(a)-EnB(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp2_Temp(2) = Temp1  
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do k = 1,NoccB
      Temp2 = 0.0d0
      do l = 1,NoccA
        Temp3 = 0.0d0
        do a = NoccB+1,N
          Temp4 = 0.0d0
          do b = NoccA+1,N        
             Temp4 = Temp4 + (g_MoAB(EI(a,k),EI(b,l)))**2/(EnB(k)+EnA(l)-EnB(a)-EnA(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp2_Temp(3) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do k = 1,NoccA
      Temp2 = 0.0d0
      do l = 1,NoccB
        Temp3 = 0.0d0
        do a = NoccB+1,N
          Temp4 = 0.0d0
          do b = NoccA+1,N        
             Temp4 = Temp4 + g_MoAB(EI(a,l),EI(b,k))*g_MoAB(EI(a,l),EI(b,k))/(EnA(k)+EnB(l)-EnB(a)-EnA(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp2_Temp(4) = Temp1  
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do k = 1,NoccA
      Temp2 = 0.0d0
      do l = 1,NoccA
        Temp3 = 0.0d0
        do a = NoccA+1,N
          Temp4 = 0.0d0
          do b = NoccA+1,N        
             Temp4 = Temp4 + (g_MoA(EI(b,l),EI(a,k))-g_MoA(EI(a,l),EI(b,k)))**2/(EnA(k)+EnA(l)-EnA(a)-EnA(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp2_Temp(5) = Temp1  
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do k = 1,NoccA
      Temp2 = 0.0d0
      do l = 1,NoccB
        Temp3 = 0.0d0
        do a = NoccA+1,N
          Temp4 = 0.0d0
          do b = NoccB+1,N        
             Temp4 = Temp4 + (g_MoAB(EI(b,l),EI(a,k)))**2/(EnA(k)+EnB(l)-EnA(a)-EnB(b))
          end do
          Temp3 = Temp3 + Temp4
        end do
        Temp2 = Temp2 + Temp3
      end do
      Temp1 = Temp1 + Temp2
    end do
    
    mp2_Temp(6) = Temp1 
    !-------------------------------------------------------------------
    mp_2 = 0.25*(mp2_Temp(1)+mp2_Temp(2)+mp2_Temp(3)+mp2_Temp(4)+mp2_Temp(5)+mp2_Temp(6))  
    
end subroutine MP2_Open
!=======================================================================
subroutine MP3_Open(N,NoccA,NoccB,EnA,EnB,g_MoA,g_MoB,g_MoAB,EI,mp_3)
    !------------- Open Shell,Multiplicity > 1 --------------!
    implicit none
    integer N,NoccA,NoccB,i,EI(N,N)
    integer a,b,c,d,r,s,t,u
    real(kind=8) mp_3,mp3_Temp(10)
    real(kind=8) EnA(N),EnB(N),g_MoAB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) Temp1,Temp2,Temp3,Temp4,Temp5,Temp6
    real(kind=8) R1,R2
    real(kind=8) e3_I,e3_II,e3_III

    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do d = 1,NoccA
        Temp5 = 0.0d0
        do r = NoccB+1,N
         Temp6 = 0.0d0
         do s = NoccA+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnB(c)+EnA(d)-EnB(r)-EnA(s))
            Temp6 = Temp6 + g_MoAB(EI(r,a),EI(s,b))*g_MoAB(EI(a,c),EI(b,d))*g_MoAB(EI(c,r),EI(d,s))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(1) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do d = 1,NoccA
        Temp5 = 0.0d0
        do r = NoccB+1,N
         Temp6 = 0.0d0
         do s = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnB(c)+EnA(d)-EnB(r)-EnA(s))
            Temp6 = Temp6 + g_MoAB(EI(r,b),EI(s,a))*g_MoAB(EI(b,c),EI(a,d))*g_MoAB(EI(c,r),EI(d,s))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(2) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do d = 1,NoccB
        Temp5 = 0.0d0
        do r = NoccB+1,N
         Temp6 = 0.0d0
         do s = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnB(b)-EnB(r)-EnB(s))
            R2 = 1.0d0/(EnB(c)+EnB(d)-EnB(r)-EnB(s))
            Temp6 = Temp6 + R1*R2*((g_MoB(EI(s,b),EI(r,a))-g_MoB(EI(r,b),EI(s,a)))* &
                    (g_MoB(EI(b,d),EI(a,c))-g_MoB(EI(a,d),EI(b,c)))* &
                    (g_MoB(EI(d,s),EI(c,r))-g_MoB(EI(c,s),EI(d,r))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(3) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do d = 1,NoccB
        Temp5 = 0.0d0
        do r = NoccB+1,N
         Temp6 = 0.0d0
         do s = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnA(c)+EnB(d)-EnB(r)-EnA(s))
            Temp6 = Temp6 + g_MoAB(EI(r,b),EI(s,a))*g_MoAB(EI(b,d),EI(a,c))*g_MoAB(EI(d,r),EI(c,s))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(4) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do d = 1,NoccB
        Temp5 = 0.0d0
        do r = NoccB+1,N
         Temp6 = 0.0d0
         do s = NoccA+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnA(c)+EnB(d)-EnB(r)-EnA(s))
            Temp6 = Temp6 + g_MoAB(EI(r,a),EI(s,b))*g_MoAB(EI(a,d),EI(b,c))*g_MoAB(EI(d,r),EI(c,s))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(5) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do d = 1,NoccB
        Temp5 = 0.0d0
        do r = NoccA+1,N
         Temp6 = 0.0d0
         do s = NoccB+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnA(c)+EnB(d)-EnA(r)-EnB(s))
            Temp6 = Temp6 + g_MoAB(EI(s,b),EI(r,a))*g_MoAB(EI(b,d),EI(a,c))*g_MoAB(EI(d,s),EI(c,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(6) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do d = 1,NoccB
        Temp5 = 0.0d0
        do r = NoccA+1,N
         Temp6 = 0.0d0
         do s = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnA(c)+EnB(d)-EnA(r)-EnB(s))
            Temp6 = Temp6 + g_MoAB(EI(s,a),EI(r,b))*g_MoAB(EI(a,d),EI(b,c))*g_MoAB(EI(d,s),EI(c,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(7) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do d = 1,NoccA
        Temp5 = 0.0d0
        do r = NoccA+1,N
         Temp6 = 0.0d0
         do s = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnB(c)+EnA(d)-EnA(r)-EnB(s))
            Temp6 = Temp6 + g_MoAB(EI(s,a),EI(r,b))*g_MoAB(EI(a,c),EI(b,d))*g_MoAB(EI(c,s),EI(d,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(8) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do d = 1,NoccA
        Temp5 = 0.0d0
        do r = NoccA+1,N
         Temp6 = 0.0d0
         do s = NoccB+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnB(c)+EnA(d)-EnA(r)-EnB(s))
            Temp6 = Temp6 + g_MoAB(EI(s,b),EI(r,a))*g_MoAB(EI(b,c),EI(a,d))*g_MoAB(EI(c,s),EI(d,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(9) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do d = 1,NoccA
        Temp5 = 0.0d0
        do r = NoccA+1,N
         Temp6 = 0.0d0
         do s = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnA(b)-EnA(r)-EnA(s))
            R2 = 1.0d0/(EnA(c)+EnA(d)-EnA(r)-EnA(s))
            Temp6 = Temp6 + R1*R2*((g_MoA(EI(s,b),EI(r,a))-g_MoA(EI(r,b),EI(s,a)))* &
                    (g_MoA(EI(b,d),EI(a,c))-g_MoA(EI(a,d),EI(b,c)))*(g_MoA(EI(d,s),EI(c,r))-g_MoA(EI(c,s),EI(d,r))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(10) = Temp1
    !-------------------------------------------------------------------
    e3_I = 0.0d0
    do i =1,10
        e3_I = e3_I + mp3_Temp(i)
    end do
    e3_I = 0.125d0*e3_I
   !====================================================================
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do r = NoccB+1,N
       Temp4 = 0.0d0
       do s = NoccA+1,N
        Temp5 = 0.0d0
        do t = NoccB+1,N
         Temp6 = 0.0d0
         do u = NoccA+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnB(a)+EnA(b)-EnB(t)-EnA(u))
            Temp6 = Temp6 + g_MoAB(EI(r,a),EI(s,b))*g_MoAB(EI(t,r),EI(u,s))*g_MoAB(EI(a,t),EI(b,u))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(1) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do r = NoccA+1,N
       Temp4 = 0.0d0
       do s = NoccB+1,N
        Temp5 = 0.0d0
        do t = NoccB+1,N
         Temp6 = 0.0d0
         do u = NoccA+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnB(a)+EnA(b)-EnB(t)-EnA(u))
            Temp6 = Temp6 + g_MoAB(EI(s,a),EI(r,b))*g_MoAB(EI(t,s),EI(u,r))*g_MoAB(EI(a,t),EI(b,u))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(2) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do r = NoccB+1,N
       Temp4 = 0.0d0
       do s = NoccB+1,N
        Temp5 = 0.0d0
        do t = NoccB+1,N
         Temp6 = 0.0d0
         do u = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnB(b)-EnB(r)-EnB(s))
            R2 = 1.0d0/(EnB(a)+EnB(b)-EnB(t)-EnB(u))
            Temp6 = Temp6 + R1*R2*((g_MoB(EI(s,b),EI(r,a))-g_MoB(EI(r,b),EI(s,a)))* &
                    (g_MoB(EI(u,s),EI(t,r))-g_MoB(EI(t,s),EI(u,r)))* &
                    (g_MoB(EI(b,u),EI(a,t))-g_MoB(EI(a,u),EI(b,t))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(3) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do r = NoccB+1,N
       Temp4 = 0.0d0
       do s = NoccA+1,N
        Temp5 = 0.0d0
        do t = NoccB+1,N
         Temp6 = 0.0d0
         do u = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnA(a)+EnB(b)-EnB(t)-EnA(u))
            Temp6 = Temp6 + g_MoAB(EI(r,b),EI(s,a))*g_MoAB(EI(t,r),EI(u,s))*g_MoAB(EI(b,t),EI(a,u))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(4) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do r = NoccA+1,N
       Temp4 = 0.0d0
       do s = NoccB+1,N
        Temp5 = 0.0d0
        do t = NoccB+1,N
         Temp6 = 0.0d0
         do u = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnA(a)+EnB(b)-EnB(t)-EnA(u))
            Temp6 = Temp6 + g_MoAB(EI(s,b),EI(r,a))*g_MoAB(EI(t,s),EI(u,r))*g_MoAB(EI(b,t),EI(a,u))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(5) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do r = NoccB+1,N
       Temp4 = 0.0d0
       do s = NoccA+1,N
        Temp5 = 0.0d0
        do t = NoccA+1,N
         Temp6 = 0.0d0
         do u = NoccB+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnA(a)+EnB(b)-EnA(t)-EnB(u))
            Temp6 = Temp6 + g_MoAB(EI(r,b),EI(s,a))*g_MoAB(EI(u,r),EI(t,s))*g_MoAB(EI(b,u),EI(a,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(6) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do r = NoccA+1,N
       Temp4 = 0.0d0
       do s = NoccB+1,N
        Temp5 = 0.0d0
        do t = NoccA+1,N
         Temp6 = 0.0d0
         do u = NoccB+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnA(a)+EnB(b)-EnA(t)-EnB(u))
            Temp6 = Temp6 + g_MoAB(EI(s,b),EI(r,a))*g_MoAB(EI(u,s),EI(t,r))*g_MoAB(EI(b,u),EI(a,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(7) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do r = NoccB+1,N
       Temp4 = 0.0d0
       do s = NoccA+1,N
        Temp5 = 0.0d0
        do t = NoccA+1,N
         Temp6 = 0.0d0
         do u = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnB(a)+EnA(b)-EnA(t)-EnB(u))
            Temp6 = Temp6 + g_MoAB(EI(r,a),EI(s,b))*g_MoAB(EI(u,r),EI(t,s))*g_MoAB(EI(a,u),EI(b,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(8) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do r = NoccA+1,N
       Temp4 = 0.0d0
       do s = NoccB+1,N
        Temp5 = 0.0d0
        do t = NoccA+1,N
         Temp6 = 0.0d0
         do u = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnB(a)+EnA(b)-EnA(t)-EnB(u))
            Temp6 = Temp6 + g_MoAB(EI(s,a),EI(r,b))*g_MoAB(EI(u,s),EI(t,r))*g_MoAB(EI(a,u),EI(b,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(9) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do r = NoccA+1,N
       Temp4 = 0.0d0
       do s = NoccA+1,N
        Temp5 = 0.0d0
        do t = NoccA+1,N
         Temp6 = 0.0d0
         do u = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnA(b)-EnA(r)-EnA(s))
            R2 = 1.0d0/(EnA(a)+EnA(b)-EnA(t)-EnA(u))
            Temp6 = Temp6 + R1*R2*((g_MoA(EI(s,b),EI(r,a))-g_MoA(EI(r,b),EI(s,a)))* &
                    (g_MoA(EI(u,s),EI(t,r))-g_MoA(EI(t,s),EI(u,r)))*(g_MoA(EI(b,u),EI(a,t))-g_MoA(EI(a,u),EI(b,t))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(10) = Temp1
    !-------------------------------------------------------------------
    e3_II = 0.0d0
    do i =1,10
        e3_II = e3_II + mp3_Temp(i)
    end do
    e3_II = 0.125d0*e3_II
   !====================================================================
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do r = NoccB+1,N
        Temp5 = 0.0d0
        do s = NoccB+1,N
         Temp6 = 0.0d0
         do t = NoccA+1,N
            R1 = 1.0d0/(EnB(a)+EnB(b)-EnB(r)-EnB(s))
            R2 = 1.0d0/(EnB(a)+EnA(c)-EnB(r)-EnA(t))
            Temp6 = Temp6 + (g_MoB(EI(s,b),EI(r,a))-g_MoB(EI(r,b),EI(s,a)))*g_MoAB(EI(b,s),EI(t,c))* &
                    g_MoAB(EI(a,r),EI(c,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(1) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do r = NoccB+1,N
        Temp5 = 0.0d0
        do s = NoccA+1,N
         Temp6 = 0.0d0
         do t = NoccA+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnB(a)+EnA(c)-EnB(r)-EnA(t))
            Temp6 = Temp6 + (g_MoA(EI(b,s),EI(t,c))-g_MoA(EI(t,s),EI(b,c)))*g_MoAB(EI(r,a),EI(s,b))* &
                    g_MoAB(EI(a,r),EI(c,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(2) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do r = NoccB+1,N
        Temp5 = 0.0d0
        do s = NoccB+1,N
         Temp6 = 0.0d0
         do t = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnB(b)-EnB(r)-EnB(s))
            R2 = 1.0d0/(EnB(a)+EnB(c)-EnB(r)-EnB(t))
            Temp6 = Temp6 + R1*R2*((g_MoB(EI(s,b),EI(r,a))-g_MoB(EI(r,b),EI(s,a)))* &
                    (g_MoB(EI(b,s),EI(t,c))-g_MoB(EI(t,s),EI(b,c)))* &
                    (g_MoB(EI(c,t),EI(a,r))-g_MoB(EI(a,t),EI(c,r))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(3) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do r = NoccB+1,N
        Temp5 = 0.0d0
        do s = NoccA+1,N
         Temp6 = 0.0d0
         do t = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnB(a)+EnB(c)-EnB(r)-EnB(t))
            Temp6 = Temp6 + (g_MoB(EI(c,t),EI(a,r))-g_MoB(EI(a,t),EI(c,r)))*g_MoAB(EI(r,a),EI(s,b))* &
                    g_MoAB(EI(t,c),EI(b,s))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(4) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do r = NoccB+1,N
        Temp5 = 0.0d0
        do s = NoccA+1,N
         Temp6 = 0.0d0
         do t = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnB(r)-EnA(s))
            R2 = 1.0d0/(EnA(a)+EnB(c)-EnB(r)-EnA(t))
            Temp6 = Temp6 + g_MoAB(EI(r,b),EI(s,a))*g_MoAB(EI(b,c),EI(t,s))*g_MoAB(EI(c,r),EI(a,t))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(5) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do r = NoccA+1,N
        Temp5 = 0.0d0
        do s = NoccB+1,N
         Temp6 = 0.0d0
         do t = NoccB+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnA(a)+EnB(c)-EnA(r)-EnB(t))
            Temp6 = Temp6 + (g_MoB(EI(b,s),EI(t,c))-g_MoB(EI(t,s),EI(b,c)))*g_MoAB(EI(s,b),EI(r,a))* &
                    g_MoAB(EI(c,t),EI(a,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(6) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccB
       Temp4 = 0.0d0
       do r = NoccA+1,N
        Temp5 = 0.0d0
        do s = NoccA+1,N
         Temp6 = 0.0d0
         do t = NoccB+1,N
            R1 = 1.0d0/(EnA(a)+EnA(b)-EnA(r)-EnA(s))
            R2 = 1.0d0/(EnA(a)+EnB(c)-EnA(r)-EnB(t))
            Temp6 = Temp6 + (g_MoA(EI(s,b),EI(r,a))-g_MoA(EI(r,b),EI(s,a)))*g_MoAB(EI(t,c),EI(b,s))* &
                    g_MoAB(EI(c,t),EI(a,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(7) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccB
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do r = NoccA+1,N
        Temp5 = 0.0d0
        do s = NoccB+1,N
         Temp6 = 0.0d0
         do t = NoccB+1,N
            R1 = 1.0d0/(EnB(a)+EnA(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnB(a)+EnA(c)-EnA(r)-EnB(t))
            Temp6 = Temp6 + g_MoAB(EI(s,a),EI(r,b))*g_MoAB(EI(t,s),EI(b,c))*g_MoAB(EI(a,t),EI(c,r))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(8) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccB
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do r = NoccA+1,N
        Temp5 = 0.0d0
        do s = NoccB+1,N
         Temp6 = 0.0d0
         do t = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnB(b)-EnA(r)-EnB(s))
            R2 = 1.0d0/(EnA(a)+EnA(c)-EnA(r)-EnA(t))
            Temp6 = Temp6 + (g_MoA(EI(c,t),EI(a,r))-g_MoA(EI(a,t),EI(c,r)))*g_MoAB(EI(s,b),EI(r,a))* &
                    g_MoAB(EI(b,s),EI(t,c))*R1*R2
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(9) = Temp1
    !-------------------------------------------------------------------
    Temp1 = 0.0d0
    do a = 1,NoccA
     Temp2 = 0.0d0
     do b = 1,NoccA
      Temp3 = 0.0d0
      do c = 1,NoccA
       Temp4 = 0.0d0
       do r = NoccA+1,N
        Temp5 = 0.0d0
        do s = NoccA+1,N
         Temp6 = 0.0d0
         do t = NoccA+1,N
            R1 = 1.0d0/(EnA(a)+EnA(b)-EnA(r)-EnA(s))
            R2 = 1.0d0/(EnA(a)+EnA(c)-EnA(r)-EnA(t))
            Temp6 = Temp6 + R1*R2*((g_MoA(EI(s,b),EI(r,a))-g_MoA(EI(r,b),EI(s,a)))* &
                    (g_MoA(EI(b,s),EI(t,c))-g_MoA(EI(t,s),EI(b,c)))*(g_MoA(EI(c,t),EI(a,r))-g_MoA(EI(a,t),EI(c,r))))
         end do
         Temp5 = Temp5 + Temp6
        end do
        Temp4 = Temp4 + Temp5
       end do
       Temp3 = Temp3 + Temp4
      end do
      Temp2 = Temp2 + Temp3
     end do
     Temp1 = Temp1 + Temp2
    end do
    mp3_Temp(10) = Temp1
    !-------------------------------------------------------------------
    e3_III = mp3_Temp(1) + mp3_Temp(2) + mp3_Temp(3) + mp3_Temp(4) -  &
             mp3_Temp(5) + mp3_Temp(6) + mp3_Temp(7) - mp3_Temp(8) +  &
             mp3_Temp(9) + mp3_Temp(10) 
    !===================================================================
    
    mp_3 = e3_I + e3_II + e3_III
    
end subroutine MP3_Open
!=======================================================================

