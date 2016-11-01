subroutine Mol_Obt(N,Two_Electron,EI,H_Core,MoCu,f_Mo,g_Mo,YN)
    
    implicit none
    integer N,i,j,k,l,p,q,r,s,EI(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,N*(N+1)/2),H_Core(N,N),MoCu(N,N)
    real(kind=8) f_Mo(N,N),g_Mo(N*(N+1)/2,N*(N+1)/2)
    real(kind=8),allocatable::A_T(:,:,:,:),B_T(:,:,:,:),C_T(:,:,:,:)
    real(kind=8) Temp1,Temp2
    character(len=1) YN
    
    
    if( YN == "Y" .or. YN == "W" ) then
      do i = 1,N
        do j = i,N
          Temp1 = 0.0d0
          do p = 1,N
             Temp2 = 0.0d0
             do q = 1,N
                Temp2 = Temp2 + MoCu(p,i)*MoCu(q,j)*H_Core(p,q)
             end do
             Temp1 = Temp1 + Temp2
          end do
          f_Mo(j,i) = Temp1  
          f_Mo(i,j) = f_Mo(j,i)
         end do
      end do         
    end if
    
    if( YN == "W" ) then
        write(*,*) " <i|h|j> transformation complete."
    end if

    allocate(A_T(N,N,N,N))
    
    do p = 1,N
      do q = 1,N
        do r = 1,N
          do l = 1,N
            A_T(l,r,q,p) = 0.0d0
            do s = 1,N
               A_T(l,r,q,p) = A_T(l,r,q,p) + MoCu(s,l)*Two_Electron(EI(s,r),EI(q,p))
            end do
          end do
        end do
      end do
    end do
    
    allocate(B_T(N,N,N,N))
        
    do p = 1,N
      do q = 1,N
        do k = 1,N
          do l = 1,N
            B_T(l,k,q,p) = 0.0d0
            do r = 1,N
               B_T(l,k,q,p) = B_T(l,k,q,p) + MoCu(r,k)*A_T(l,r,q,p)
            end do
          end do
        end do
      end do
    end do
    
    deallocate(A_T)
    allocate(C_T(N,N,N,N)) 
       
    do p = 1,N
      do j = 1,N
        do k = 1,N
          do l = 1,N
            C_T(l,k,j,p) = 0.0d0
            do q = 1,N
               C_T(l,k,j,p) = C_T(l,k,j,p) + MoCu(q,j)*B_T(l,k,q,p)
            end do
          end do
        end do
      end do
    end do
    
    deallocate(B_T)
        
    do i = 1,N
      do j = 1,i
        do k = 1,N
          do l = 1,k
            if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
            g_Mo(EI(l,k),EI(j,i)) = 0.0d0
            do p = 1,N
               g_Mo(EI(l,k),EI(j,i)) = g_Mo(EI(l,k),EI(j,i)) + MoCu(p,i)*C_T(l,k,j,p)
            end do
            g_Mo(EI(j,i),EI(l,k)) = g_Mo(EI(l,k),EI(j,i))
          end do
        end do
      end do
    end do
    
    deallocate(C_T)
    
end subroutine Mol_Obt           
!=======================================================================
subroutine Mul_Mol_Obt(N,Two_Electron,EI,MoCuA,MoCuB,g_Mo)
    
    implicit none
    integer N,i,j,k,l,p,q,r,s,EI(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,N*(N+1)/2),MoCuA(N,N),MoCuB(N,N)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2)
    real(kind=8),allocatable::A_T(:,:,:,:),B_T(:,:,:,:),C_T(:,:,:,:)
    
    allocate(A_T(N,N,N,N))
    
    do p = 1,N
      do q = 1,N
        do r = 1,N
          do l = 1,N
            A_T(l,r,q,p) = 0.0d0
            do s = 1,N
               A_T(l,r,q,p) = A_T(l,r,q,p) + MoCuB(s,l)*Two_Electron(EI(s,r),EI(q,p))
            end do
          end do
        end do
      end do
    end do
    
    allocate(B_T(N,N,N,N))
       
    do p = 1,N
      do q = 1,N
        do k = 1,N
          do l = 1,N
            B_T(l,k,q,p) = 0.0d0
            do r = 1,N
               B_T(l,k,q,p) = B_T(l,k,q,p) + MoCuB(r,k)*A_T(l,r,q,p)
            end do
          end do
        end do
      end do
    end do
    
    deallocate(A_T)
    allocate(C_T(N,N,N,N))
           
    do p = 1,N
      do j = 1,N
        do k = 1,N
          do l = 1,N
            C_T(l,k,j,p) = 0.0d0
            do q = 1,N
               C_T(l,k,j,p) = C_T(l,k,j,p) + MoCuA(q,j)*B_T(l,k,q,p)
            end do
          end do
        end do
      end do
    end do

    deallocate(B_T)
        
    do i = 1,N
      do j = 1,i
        do k = 1,N
          do l = 1,k
            g_Mo(EI(l,k),EI(j,i)) = 0.0d0
            do p = 1,N
               g_Mo(EI(l,k),EI(j,i)) = g_Mo(EI(l,k),EI(j,i)) + MoCuA(p,i)*C_T(l,k,j,p)
            end do
          end do
        end do
      end do
    end do    
    
    deallocate(C_T)
    
end subroutine Mul_Mol_Obt           
!=======================================================================
