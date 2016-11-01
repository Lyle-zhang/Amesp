subroutine RHF_Const_G(N,Nmem,Neri,Two_Electron,EI,Density,G,t2,d2)
    !--------------------- Construct G Matrix --------------------------
    implicit none
    integer i,j,k,l,N,Nmem,Neri,EI(N,N)
    real(kind=8) G(N,N),Density(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) Temp1
    character(len=10) t2
    character(len=8) d2
    
    G = 0.0d0
       
    if ( N > Nmem ) then
       !----------------------------------------------------------------
       open(unit=1110,file="/tmp/"//d2//t2//".amesp",form="unformatted")
          
       do i = 1,N
         do j = 1,i
           do k = 1,N
             do l = 1,k
       
                if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                if(dsqrt(dabs(Two_Electron(EI(j,i),1)*Two_Electron(EI(l,k),1))) < 1.0d-8 ) cycle
        
                read(1110) Temp1
                call RHF_G_Mat(N,i,j,k,l,Temp1,G,Density)
                  
             end do
           end do
         end do
       end do
          
       close(1110)
          
       do i = 1,N
          do j = 1,i
             G(j,i) = G(i,j)
          end do
       end do
       !-------------------------------------------------------------
    else
       !-------------------------------------------------------------
       do i = 1,N
         do j = 1,i
           do k = 1,N
             do l = 1,k
       
                if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                if(dsqrt(dabs(Two_Electron(EI(j,i),EI(j,i))*Two_Electron(EI(l,k),EI(l,k)))) < 1.0d-8 ) cycle
        
                Temp1 = Two_Electron(EI(l,k),EI(j,i))
                call RHF_G_Mat(N,i,j,k,l,Temp1,G,Density)
                  
             end do
           end do
         end do
       end do

       do i = 1,N
          do j = 1,i
             G(j,i) = G(i,j)
          end do
       end do
       !-------------------------------------------------------------
    end if
 
end subroutine RHF_Const_G
!=======================================================================
subroutine UHF_Const_G(N,Nmem,Neri,Two_Electron,EI,DensityA,DensityB,Coul,ExcA,ExcB,t2,d2,CE)
    !--------------------- Construct G Matrix --------------------------
    implicit none
    integer i,j,k,l,N,Nmem,Neri,EI(N,N)
    real(kind=8) DensityA(N,N),DensityB(N,N),Coul(N,N),ExcA(N,N),ExcB(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2),Density(N,N)
    real(kind=8) Temp1
    character(len=10) t2
    character(len=8) d2
    character(len=1) CE
    
    Coul = 0.0d0 ; ExcA = 0.0d0 ; ExcB = 0.0d0
    Density = (DensityA+DensityB)/2.0d0
       
    if ( N > Nmem ) then
       !----------------------------------------------------------------
       open(unit=1110,file="/tmp/"//d2//t2//".amesp",form="unformatted")
       !----------------------------------------------------------------
       if ( CE == "C")  then         
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),1)*Two_Electron(EI(l,k),1))) < 1.0d-8 ) cycle
        
                   read(1110) Temp1
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                  
                end do
              end do
            end do
          end do          
          close(1110)         
       else !-----------------------------------------------------------
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),1)*Two_Electron(EI(l,k),1))) < 1.0d-8 ) cycle
        
                   read(1110) Temp1
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                   call Exc_Mat(N,i,j,k,l,Temp1,ExcA,DensityA)
                   call Exc_Mat(N,i,j,k,l,Temp1,ExcB,DensityB)
                  
                end do
              end do
            end do
          end do          
          close(1110) 
       end if
       !----------------------------------------------------------------          
       do i = 1,N
          do j = 1,i
             Coul(j,i) = Coul(i,j)
             ExcA(j,i) = ExcA(i,j)
             ExcB(j,i) = ExcB(i,j)
          end do
       end do
       !-------------------------------------------------------------
    else
       !-------------------------------------------------------------
       if ( CE == "C")  then         
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),EI(j,i))*Two_Electron(EI(l,k),EI(l,k)))) < 1.0d-8 ) cycle
        
                   Temp1 = Two_Electron(EI(l,k),EI(j,i))
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                  
                end do
              end do
            end do
          end do           
       else !-----------------------------------------------------------
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),EI(j,i))*Two_Electron(EI(l,k),EI(l,k)))) < 1.0d-8 ) cycle
        
                   Temp1 = Two_Electron(EI(l,k),EI(j,i))
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                   call Exc_Mat(N,i,j,k,l,Temp1,ExcA,DensityA)
                   call Exc_Mat(N,i,j,k,l,Temp1,ExcB,DensityB)
                  
                end do
              end do
            end do
          end do          
       end if
       !----------------------------------------------------------------          
       do i = 1,N
          do j = 1,i
             Coul(j,i) = Coul(i,j)
             ExcA(j,i) = ExcA(i,j)
             ExcB(j,i) = ExcB(i,j)
          end do
       end do
       !-------------------------------------------------------------
    end if
 
end subroutine UHF_Const_G
!=======================================================================
subroutine RDFT_Const_G(N,Nmem,Neri,Two_Electron,EI,Density,Coul,Exc,t2,d2,CE)
    !--------------------- Construct G Matrix --------------------------
    implicit none
    integer i,j,k,l,N,Nmem,Neri,EI(N,N)
    real(kind=8) Density(N,N),Coul(N,N),Exc(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)
    real(kind=8) Temp1
    character(len=10) t2
    character(len=8) d2
    character(len=1) CE
    
    Coul = 0.0d0 ; Exc = 0.0d0
       
    if ( N > Nmem ) then
       !----------------------------------------------------------------
       open(unit=1110,file="/tmp/"//d2//t2//".amesp",form="unformatted")
       !----------------------------------------------------------------
       if ( CE == "C")  then         
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),1)*Two_Electron(EI(l,k),1))) < 1.0d-8 ) cycle
        
                   read(1110) Temp1
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                  
                end do
              end do
            end do
          end do          
          close(1110)         
       else !-----------------------------------------------------------
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),1)*Two_Electron(EI(l,k),1))) < 1.0d-8 ) cycle
        
                   read(1110) Temp1
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                   call Exc_Mat(N,i,j,k,l,Temp1,Exc,Density)
                  
                end do
              end do
            end do
          end do          
          close(1110) 
       end if
       !----------------------------------------------------------------          
       do i = 1,N
          do j = 1,i
             Coul(j,i) = Coul(i,j)
             Exc(j,i) = Exc(i,j)
          end do
       end do
       !-------------------------------------------------------------
    else
       !-------------------------------------------------------------
       if ( CE == "C")  then         
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),EI(j,i))*Two_Electron(EI(l,k),EI(l,k)))) < 1.0d-8 ) cycle
        
                   Temp1 = Two_Electron(EI(l,k),EI(j,i))
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                  
                end do
              end do
            end do
          end do             
       else !-----------------------------------------------------------
          do i = 1,N
            do j = 1,i
              do k = 1,N
                do l = 1,k
       
                   if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
                   if(dsqrt(dabs(Two_Electron(EI(j,i),EI(j,i))*Two_Electron(EI(l,k),EI(l,k)))) < 1.0d-8 ) cycle
        
                   Temp1 = Two_Electron(EI(l,k),EI(j,i))
                   call Coul_Mat(N,i,j,k,l,Temp1,Coul,Density)
                   call Exc_Mat(N,i,j,k,l,Temp1,Exc,Density)
                  
                end do
              end do
            end do
          end do          
       end if
       !----------------------------------------------------------------          
       do i = 1,N
          do j = 1,i
             Coul(j,i) = Coul(i,j)
             Exc(j,i) = Exc(i,j)
          end do
       end do
       !-------------------------------------------------------------
    end if
 
end subroutine RDFT_Const_G
!=======================================================================
subroutine RHF_G_Mat(N,i,j,k,l,eri,G,Density)
    !-------------------- G matrix on erery eri ------------------------
    implicit none
    integer N,i,j,k,l
    real(kind=8) G(N,N),Density(N,N),eri
    !----- for (i,j|k,l) ,here i>j, i*(i-1)+2*j > k*(k-1)+2*l, k>l -----

    if (i == j .and. j == k .and. k == l) then
        G(i,i) = G(i,i) + Density(i,i)*eri
        return
    end if
    
    if (i == j .and. j == k) then
        G(i,i) = G(i,i) + 2*Density(i,l)*eri 
        G(i,l) = G(i,l) + Density(i,i)*eri
        Return 
    end if
    
    if (j == k .and. k == l) then
        G(j,j) = G(j,j) + 2*Density(i,j)*eri 
        G(i,j) = G(i,j) + Density(j,j)*eri 
        return 
    end if
    
    if (i == k .and. j == l) then
        G(i,i) = G(i,i) - Density(j,j)*eri 
        G(j,j) = G(j,j) - Density(i,i)*eri
        G(i,j) = G(i,j) + 3*Density(i,j)*eri
        return 
    end if
    
    if (i == j .and. k == l) then
        G(i,i) = G(i,i) + 2*Density(k,k)*eri
        G(k,k) = G(k,k) + 2*Density(i,i)*eri
        G(i,k) = G(i,k) - Density(i,k)*eri
        return 
    end if
          
    if (i == j) then
        G(i,i) = G(i,i) + 4*Density(k,l)*eri
        G(k,l) = G(k,l) + 2*Density(i,i)*eri
        G(i,k) = G(i,k) - Density(i,l)*eri
        G(i,l) = G(i,l) - Density(i,k)*eri
        return
    end if
           
    if (i == k) then
        G(i,j) = G(i,j) + 3*Density(i,l)*eri
        G(i,l) = G(i,l) + 3*Density(i,j)*eri
        G(j,l) = G(j,l) - Density(i,i)*eri
        G(i,i) = G(i,i) - 2*Density(j,l)*eri
        return 
    end if
         
    if (j == l) then
        G(i,j) = G(i,j) + 3*Density(k,j)*eri
        G(k,j) = G(k,j) + 3*Density(i,j)*eri
        G(j,j) = G(j,j) - 2*Density(i,k)*eri
        G(i,k) = G(i,k) - Density(j,j)*eri
        return 
    end if
         
    if (k == l  .and. j > l) then
        G(i,j) = G(i,j) + 2*Density(k,k)*eri
        G(k,k) = G(k,k) + 4*Density(i,j)*eri
        G(j,k) = G(j,k) - Density(i,k)*eri
        G(i,k) = G(i,k) - Density(j,k)*eri
        return 
    end if
    
    if (k == l) then
        G(i,j) = G(i,j) + 2*Density(k,k)*eri
        G(k,k) = G(k,k) + 4*Density(i,j)*eri
        G(k,j) = G(k,j) - Density(i,k)*eri
        G(i,k) = G(i,k) - Density(k,j)*eri
        return 
    end if
          
    if (j == k) then
        G(i,j) = G(i,j) + 3*Density(j,l)*eri
        G(j,l) = G(j,l) + 3*Density(i,j)*eri
        G(i,l) = G(i,l) - Density(j,j)*eri
        G(j,j) = G(j,j) - 2*Density(i,l)*eri
        return 
    end if     
     
    if (j > k) then
        G(i,j) = G(i,j) + 4*Density(k,l)*eri
        G(k,l) = G(k,l) + 4*Density(i,j)*eri
        G(i,l) = G(i,l) - Density(j,k)*eri
        G(j,k) = G(j,k) - Density(i,l)*eri
        G(j,l) = G(j,l) - Density(i,k)*eri
        G(i,k) = G(i,k) - Density(j,l)*eri
        return 
    end if 
        
    if (j > l) then
        G(i,j) = G(i,j) + 4*Density(k,l)*eri
        G(k,l) = G(k,l) + 4*Density(i,j)*eri
        G(i,l) = G(i,l) - Density(k,j)*eri
        G(k,j) = G(k,j) - Density(i,l)*eri
        G(j,l) = G(j,l) - Density(i,k)*eri
        G(i,k) = G(i,k) - Density(j,l)*eri
        return 
    end if
          
    if (j < l) then
        G(i,j) = G(i,j) + 4*Density(k,l)*eri
        G(k,l) = G(k,l) + 4*Density(i,j)*eri
        G(i,l) = G(i,l) - Density(k,j)*eri
        G(k,j) = G(k,j) - Density(i,l)*eri
        G(l,j) = G(l,j) - Density(i,k)*eri
        G(i,k) = G(i,k) - Density(l,j)*eri
        return 
    end if  
    
    write(*,*)
    write(*,*) " stop at G_Matrix."
    stop
      
end subroutine RHF_G_Mat
!=======================================================================
subroutine Coul_Mat(N,i,j,k,l,eri,Coul,Density)
    !------------------- Coul matrix on erery eri ----------------------
    implicit none
    integer N,i,j,k,l
    real(kind=8) Coul(N,N),Density(N,N),eri
    !----- for (i,j|k,l) ,here i>j, i*(i-1)+2*j > k*(k-1)+2*l, k>l -----

    if (i == j .and. j == k .and. k == l) then
        Coul(i,i) = Coul(i,i) + 2*Density(i,i)*eri
        return
    end if
    
    if (i == j .and. j == k) then
        Coul(i,i) = Coul(i,i) + 4*Density(i,l)*eri 
        Coul(i,l) = Coul(i,l) + 2*Density(i,i)*eri
        Return 
    end if
    
    if (j == k .and. k == l) then
        Coul(j,j) = Coul(j,j) + 4*Density(i,j)*eri 
        Coul(i,j) = Coul(i,j) + 2*Density(j,j)*eri 
        return 
    end if
    
    if (i == k .and. j == l) then
        Coul(i,j) = Coul(i,j) + 4*Density(i,j)*eri
        return 
    end if
    
    if (i == j .and. k == l) then
        Coul(i,i) = Coul(i,i) + 2*Density(k,k)*eri
        Coul(k,k) = Coul(k,k) + 2*Density(i,i)*eri
        return 
    end if
          
    if (i == j) then
        Coul(i,i) = Coul(i,i) + 4*Density(k,l)*eri
        Coul(k,l) = Coul(k,l) + 2*Density(i,i)*eri
        return
    end if
           
    if (i == k) then
        Coul(i,j) = Coul(i,j) + 4*Density(i,l)*eri
        Coul(i,l) = Coul(i,l) + 4*Density(i,j)*eri
        return 
    end if
         
    if (j == l) then
        Coul(i,j) = Coul(i,j) + 4*Density(k,j)*eri
        Coul(k,j) = Coul(k,j) + 4*Density(i,j)*eri
        return 
    end if
         
    if (k == l  .and. j > l) then
        Coul(i,j) = Coul(i,j) + 2*Density(k,k)*eri
        Coul(k,k) = Coul(k,k) + 4*Density(i,j)*eri
        return 
    end if
    
    if (k == l) then
        Coul(i,j) = Coul(i,j) + 2*Density(k,k)*eri
        Coul(k,k) = Coul(k,k) + 4*Density(i,j)*eri
        return 
    end if
          
    if (j == k) then
        Coul(i,j) = Coul(i,j) + 4*Density(j,l)*eri
        Coul(j,l) = Coul(j,l) + 4*Density(i,j)*eri
        return 
    end if     
     
    if (j > k) then
        Coul(i,j) = Coul(i,j) + 4*Density(k,l)*eri
        Coul(k,l) = Coul(k,l) + 4*Density(i,j)*eri
        return 
    end if 
        
    if (j > l) then
        Coul(i,j) = Coul(i,j) + 4*Density(k,l)*eri
        Coul(k,l) = Coul(k,l) + 4*Density(i,j)*eri
        return 
    end if
          
    if (j < l) then
        Coul(i,j) = Coul(i,j) + 4*Density(k,l)*eri
        Coul(k,l) = Coul(k,l) + 4*Density(i,j)*eri
        return 
    end if  
    
    write(*,*)
    write(*,*) " stop at Coul_Matrix."
    stop
      
end subroutine Coul_Mat
!=======================================================================
subroutine Exc_Mat(N,i,j,k,l,eri,Exc,Density)
    !-------------------- Exc matrix on erery eri ----------------------
    implicit none
    integer N,i,j,k,l
    real(kind=8) Exc(N,N),Density(N,N),eri
    !----- for (i,j|k,l) ,here i>j, i*(i-1)+2*j > k*(k-1)+2*l, k>l -----

    if (i == j .and. j == k .and. k == l) then
        Exc(i,i) = Exc(i,i) - Density(i,i)*eri
        return
    end if
    
    if (i == j .and. j == k) then
        Exc(i,i) = Exc(i,i) - 2*Density(i,l)*eri 
        Exc(i,l) = Exc(i,l) - Density(i,i)*eri
        Return 
    end if
    
    if (j == k .and. k == l) then
        Exc(j,j) = Exc(j,j) - 2*Density(i,j)*eri 
        Exc(i,j) = Exc(i,j) - Density(j,j)*eri 
        return 
    end if
    
    if (i == k .and. j == l) then
        Exc(i,i) = Exc(i,i) - Density(j,j)*eri 
        Exc(j,j) = Exc(j,j) - Density(i,i)*eri
        Exc(i,j) = Exc(i,j) - Density(i,j)*eri
        return 
    end if
    
    if (i == j .and. k == l) then
        Exc(i,k) = Exc(i,k) - Density(i,k)*eri
        return 
    end if
          
    if (i == j) then
        Exc(i,k) = Exc(i,k) - Density(i,l)*eri
        Exc(i,l) = Exc(i,l) - Density(i,k)*eri
        return
    end if
           
    if (i == k) then
        Exc(i,j) = Exc(i,j) - Density(i,l)*eri
        Exc(i,l) = Exc(i,l) - Density(i,j)*eri
        Exc(j,l) = Exc(j,l) - Density(i,i)*eri
        Exc(i,i) = Exc(i,i) - 2*Density(j,l)*eri
        return 
    end if
         
    if (j == l) then
        Exc(i,j) = Exc(i,j) - Density(k,j)*eri
        Exc(k,j) = Exc(k,j) - Density(i,j)*eri
        Exc(j,j) = Exc(j,j) - 2*Density(i,k)*eri
        Exc(i,k) = Exc(i,k) - Density(j,j)*eri
        return 
    end if
         
    if (k == l  .and. j > l) then
        Exc(j,k) = Exc(j,k) - Density(i,k)*eri
        Exc(i,k) = Exc(i,k) - Density(j,k)*eri
        return 
    end if
    
    if (k == l) then
        Exc(k,j) = Exc(k,j) - Density(i,k)*eri
        Exc(i,k) = Exc(i,k) - Density(k,j)*eri
        return 
    end if
          
    if (j == k) then
        Exc(i,j) = Exc(i,j) - Density(j,l)*eri
        Exc(j,l) = Exc(j,l) - Density(i,j)*eri
        Exc(i,l) = Exc(i,l) - Density(j,j)*eri
        Exc(j,j) = Exc(j,j) - 2*Density(i,l)*eri
        return 
    end if     
     
    if (j > k) then
        Exc(i,l) = Exc(i,l) - Density(j,k)*eri
        Exc(j,k) = Exc(j,k) - Density(i,l)*eri
        Exc(j,l) = Exc(j,l) - Density(i,k)*eri
        Exc(i,k) = Exc(i,k) - Density(j,l)*eri
        return 
    end if 
        
    if (j > l) then
        Exc(i,l) = Exc(i,l) - Density(k,j)*eri
        Exc(k,j) = Exc(k,j) - Density(i,l)*eri
        Exc(j,l) = Exc(j,l) - Density(i,k)*eri
        Exc(i,k) = Exc(i,k) - Density(j,l)*eri
        return 
    end if
          
    if (j < l) then
        Exc(i,l) = Exc(i,l) - Density(k,j)*eri
        Exc(k,j) = Exc(k,j) - Density(i,l)*eri
        Exc(l,j) = Exc(l,j) - Density(i,k)*eri
        Exc(i,k) = Exc(i,k) - Density(l,j)*eri
        return 
    end if  
    
    write(*,*)
    write(*,*) " stop at Exc_Matrix."
    stop
      
end subroutine Exc_Mat
!=======================================================================
subroutine Nucleus_Energy(M,Rn,MAtom,E_Nuc)           ! Energy between Nucleus

    implicit none
    integer M,i,j,MAtom(M)
    real(kind=8) Rn(3,M),Temp1,Temp2,E_Nuc
    
    if( M == 1) then
         E_Nuc = 0.0d0
         return
    end if
    
    Temp1 = 0.0d0
    do i = 1,M
        Temp2 = 0.0d0
        do j = i+1,M
            Temp2 = Temp2 + MAtom(i)*MAtom(j)/dsqrt((Rn(1,i)-Rn(1,j))**2+(Rn(2,i)- &
                    Rn(2,j))**2+(Rn(3,i)-Rn(3,j))**2)
        end do
        Temp1 = Temp1 + Temp2
    end do
    E_Nuc = Temp1 
    
    if( M == 2) then
         E_Nuc = MAtom(1)*MAtom(2)/dsqrt((Rn(1,1)-Rn(1,2))**2+(Rn(2,1)- &
                 Rn(2,2))**2+(Rn(3,1)-Rn(3,2))**2)  
    end if           
                                               
end subroutine Nucleus_Energy
!==============================DIIS=====================================
subroutine Error_Mat(N,Fn,S,Pn_1,Er,emax)

    implicit none
    integer i,j,N
    real(kind=8) emax
    real(kind=8) Fn(N,N),S(N,N),Pn_1(N,N),Er(N,N)

    ! Error matrix, E=FDS-SDF   
    Er = matmul(matmul(Fn,Pn_1),S) - matmul(matmul(S,Pn_1),Fn)
    
    ! Find the max element of Error matrix.
    emax = 0.0d0
    do i = 1,N
       do j = 1,N
          emax = max(emax,dabs(Er(i,j)))
       end do
    end do
       
end subroutine Error_Mat
!==============================DIIS=====================================
subroutine C_Coff(ii,N,Er,CCoff)

    implicit none
    integer i,j,k,ii,N
    real(kind=8) Er(N,N,ii),CCoff(ii+1),B(ii+1,ii+1),Bb(ii+1)
    real(kind=8) MTemp(N,N)
    
    do i = 1,ii
       do j = 1,ii
          MTemp = matmul(Er(:,:,i),Er(:,:,j))
          B(i,j) = 0.0d0
          do k = 1,N
              B(i,j) = B(i,j) + MTemp(k,k)
          end do
       end do
    end do
    
    B(ii+1,1:ii) = -1.0d0
    B(1:ii,ii+1) = -1.0d0
    B(ii+1,ii+1) = 0.0d0
    
    Bb = 0.0d0 ; Bb(ii+1) = -1.0d0
    
    call LinearSolve(B,Bb,CCoff,ii+1)
    
end subroutine C_Coff
!==============================DIIS=====================================
subroutine UC_Coff(ii,N,ErA,ErB,CCoff)

    implicit none
    integer i,j,k,ii,N
    real(kind=8) ErA(N,N,ii),ErB(N,N,ii),CCoff(ii+1),B(ii+1,ii+1),Bb(ii+1)
    real(kind=8) MTemp(N,N)
    
    do i = 1,ii
       do j = 1,ii
          MTemp = matmul(ErA(:,:,i),ErB(:,:,j))
          B(i,j) = 0.0d0
          do k = 1,N
              B(i,j) = B(i,j) + MTemp(k,k)
          end do
       end do
    end do
    
    B(ii+1,1:ii) = -1.0d0
    B(1:ii,ii+1) = -1.0d0
    B(ii+1,ii+1) = 0.0d0
    
    Bb = 0.0d0 ; Bb(ii+1) = -1.0d0
    
    call LinearSolve(B,Bb,CCoff,ii+1)
    
end subroutine UC_Coff
!=======================================================================
subroutine U_Spin_S2(N,NoccA,NoccB,Overlap,MoCuA,MoCuB,S2)
    !------------------- <S**2> --------------------!
    implicit none
    integer i,j,k,l,N,NoccA,NoccB
    real(kind=8) Overlap(N,N),MoCuA(N,N),MoCuB(N,N),S2
    real(kind=8) Temp1,S2_Temp
    
    S2_Temp = (NoccA-NoccB)/2.0d0
    S2 = S2_Temp*(S2_Temp + 1.0d0) + NoccB
    
    do i = 1,NoccA
        do j = 1,NoccB
           Temp1 = 0.0d0
           do k = 1,N
               do l = 1,N
                   Temp1 = Temp1 + MoCuA(k,i)*MoCuB(l,j)*Overlap(k,l)
               end do
           end do
           S2 = S2 - Temp1**2
        end do
    end do

end subroutine U_Spin_S2
!=======================================================================
subroutine Main_Calculate_S(Sqrt_S,F,C,E,N)
!--------Solve FC=SCE by Symmetric Orthogonalization--------!

    implicit none
    integer N,i,j
    real(kind=8) E(N)
    real(kind=8) F(N,N),C(N,N),Temp(N,N)
    real(kind=8) Sqrt_S(N,N),C_Temp(N,N),F_Temp(N,N)
    
    Temp = matmul(Sqrt_S,F)
    F_Temp = matmul(Temp,Sqrt_S) ! Calculate F' = S^(-1/2)FS^(-1/2)
    
    do i = 1,N
       do j = 1,N
           if( dabs(F_Temp(i,j)) < 1d-14 ) F_Temp(i,j) = dabs(F_Temp(i,j))
       end do
    end do
    !-------------------------------------------------------------------
    C_Temp = 0.0d0
    call SY_Eigen_System(F_Temp,E,C_Temp,N)  ! Solve F'C'=EC'

    do i = 1,N
       do j = 1,N
           if( dabs(C_Temp(i,j)) < 1d-14 ) C_Temp(i,j) = dabs(C_Temp(i,j))
       end do
    end do
    !-------------------------------------------------------------------
    C = matmul(Sqrt_S,C_Temp)                ! Calculate C = S^(-1/2)*C'
    
    do i = 1,N
       do j = 1,N
           if( dabs(C(i,j)) < 1d-14 ) C(i,j) = dabs(C(i,j))
       end do
    end do
      
end subroutine Main_Calculate_S
!=======================================================================
subroutine Main_Calculate_C(S,F,C,E,N)
!--------Solve FC=SCE by Canonical Orthogonalization--------!

    implicit none
    integer N,i
    real(kind=8) E(N),Temp(N,N)
    real(kind=8) S(N,N),F(N,N),C(N,N)
    real(kind=8) X(N,N),C_Temp(N,N),F_Temp(N,N)
    
    call SY_Eigen_System(S,E,X,N)              ! Diagonalizing S
    
    do i = 1,N
         X(:,i) = X(:,i)/dsqrt(E(i))
    end do                                     ! Calaulate X = Us^(-1/2)
    
    Temp = matmul(transpose(X),F)
    F_Temp = matmul(Temp,X)  ! Calculate F' = X'FX
    call SY_Eigen_System(F_Temp,E,C_Temp,N)    ! Solve F'C'=EC'
    C = matmul(X,C_Temp)                       ! Calculate C = X*C'
      
end subroutine Main_Calculate_C
!=======================================================================
subroutine Mat_Sqrt(A,An,N)  
    ! An = A^(-1/2)
    
    implicit none
    integer N,i,j,k
    real(kind=8) lambda(N),Z(N,N),A(N,N)
    real(kind=8) ZM(N,N),An(N,N)
   
    call SY_Eigen_System(A,lambda,Z,N)

    do i = 1 ,N
         do j = 1,N
              ZM(i,j) = Z(i,j)/dsqrt(lambda(j))
         end do
    end do
   
    do i = 1,N
        do j = 1,N
            An(i,j) = 0.0d0
            do k = 1,N
                An(i,j) = An(i,j) + ZM(i,k)*Z(j,k)
            end do
            if( dabs(An(i,j)) < 1d-14 ) An(i,j) = dabs(An(i,j))
        end do
    end do
   
end subroutine Mat_Sqrt
!=======================================================================
subroutine Mat_Sqrt1(A,An,N)  
    ! An = A^(1/2)
    
    implicit none
    integer N,i,j,k
    real(kind=8) lambda(N),Z(N,N),A(N,N)
    real(kind=8) ZM(N,N),An(N,N)
   
    call SY_Eigen_System(A,lambda,Z,N)

    do i = 1 ,N
         do j = 1,N
              ZM(i,j) = Z(i,j)*dsqrt(lambda(j))
         end do
    end do
   
    do i = 1,N
        do j = 1,N
            An(i,j) = 0.0d0
            do k = 1,N
                An(i,j) = An(i,j) + ZM(i,k)*Z(j,k)
            end do
            if( dabs(An(i,j)) < 1d-14 ) An(i,j) = dabs(An(i,j))
        end do
    end do
   
end subroutine Mat_Sqrt1
!=======================================================================
subroutine SY_Eigen_System(A,lambda,Z,N)
    !Solve A*Z = lambda*Z ,A is real symmetric matrix.
    
    implicit none 
    integer N,info
    real(kind=8) A(N,N),Z(N,N),lambda(N),work(3*N-1)
    
    Z = A
    
    call dsyev('V','L', N, Z, N, lambda, work, 3*N-1, info)
    if(info /= 0) write(*,*)"segmentation fault !"
    
end subroutine SY_Eigen_System
!=======================================================================
subroutine LinearSolve(A,b,x,N)
    !Solve Ax=b
    implicit none
    integer N,info,lda,ldb
    real(kind=8) A(N,N),b(N),x(N),ipiv(N)
   
    lda = N
    ldb = N

    call dgesv (N, 1, A, lda, ipiv, b, ldb, info)
    if(info /= 0) write(*,*) "segmentation fault !"
   
    x = b

 end subroutine LinearSolve
!=======================================================================
subroutine Inverse(A,invA,N)
    !A*invA=E
    
    implicit none
    integer N,lda,ipiv(N),info,lwork
    real(kind=8) A(N,N),invA(N,N),work(N)
  
    lwork = N
    lda = N
    invA = A
   
    call dgetrf(N, N, invA, lda, ipiv, info)
    if(info/=0) write(*,*) 'Error occured in zgetrf!'
    call dgetri(N, invA, lda, ipiv, work, lwork, info)
    if(info/=0) write(*,*) 'Error occured in zgetri!'

end subroutine Inverse
!=======================================================================
