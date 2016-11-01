subroutine Match(N,str1,str2,strx,stry,strm,strn,ndif,ip)

    !  str1: State1  ;  str2: State2.
    !  strx: Diffrerent elements of State1  ;  stry: Diffrerent elements of State2.
    !  strm: Position of the different elements of State1.
    !  strn: Position of the different elements of State2.
    !  ip  : Sign of the two states after match.
    !  ndif: Number of different elements.
    
    implicit none
    integer N
    integer i,j,k,ndif,j1,ip
    integer str1(N),str2(N),strx(N),stry(N),strm(N),strn(N)

    ndif = 0 ; j1 = 1
    strx = 0 ; stry = 0 ; strm = 0 ; strn = 0
    
    do i = 1,N
       k = 0
       do j = j1,N

          if( str1(i) == str2(j) )then
              k = k + 1
              j1 = j
              exit
          end if 

       end do
       if( k == 0 ) then
           ndif = ndif + 1
           strx(ndif) = str1(i)
           strm(ndif) = i
       end if
       
    end do
    !-------------------------------------------------------------------

    ndif = 0 ; j1 = 1;
    
    do i = 1,N
       k = 0
       do j = j1,N

          if( str2(i) == str1(j) )then
              k = k + 1
              j1 = j
              exit
          end if 

       end do
       if( k == 0 ) then
           ndif = ndif + 1
           stry(ndif) = str2(i)
           strn(ndif) = i
       end if
       
    end do
    !-------------------------------------------------------------------
    j = 0 ; ip = 1
    do i = ndif,1,-1
         ip = (-1)**(N - j - strm(i))*ip
         j = j + 1
    end do
    
    j = 0
    do i = ndif,1,-1
         ip = (-1)**(N - j - strn(i))*ip
         j = j + 1
    end do
    
end subroutine Match
!=======================================================================
subroutine Bubble_Sort(A,N,ip)

    implicit none  
    integer  N,i,j
    integer  A(N),ip,temp
    
    ip = 1
    
    do i=N-1,1,-1 
       do j=1,i
          if ( A(j) > A(j+1) ) then
             ip = ip*(-1)
             temp = A(j)
             A(j) = A(j+1)
             A(j+1) = temp
          end if
       end do
    end do

end subroutine Bubble_Sort
!=======================================================================
subroutine Dig_CI_Ham(H,E,N,E_Nuc)
    ! diabonalization the ci Hamiltonian matrix.
    implicit none
    integer N,ii,i,j
    real(kind=8) E,E_Nuc,m,Temp1,ET
    real(kind=8) H(N,N),u(N),v(N)
    
    u = 0.0d0 ; u(1) = 1.0d0 ; ET = 100.0d0
    write(*,*)
    write(100,*)
    write(*,*) "************* Configuration Interaction *************"
    write(100,*) "************* Configuration Interaction *************"
    write(*,*)
    write(100,*)
    write(*,"(A7,A17,A23)") "cycle","E(CI)","Delta_E(CI)" 
    write(100,"(A7,A17,A23)") "cycle","E(CI)","Delta_E(CI)"    
    
    do ii = 1,10000*N

       Temp1 = 0.0d0
       do i = 1,N
          Temp1 = Temp1 + u(i)**2
       end do
    
       m = dsqrt(Temp1)    
       v = u/m
    
    
       do i = 1,N
          u(i) = 0.0d0
          do j = 1,N
             u(i) = u(i) + H(i,j)*v(j)
          end do
       end do
    
       E = 0.0d0
       do i = 1,N
          E = E + u(i)*v(i)
       end do
       
       write(*,"(I6,F20.7,ES22.7)") ii,E + E_Nuc,dabs(E - ET)
       
       if(dabs(Et - E) < 1.0d-9 .or. mod(ii,50) == 0 ) then
           write(100,"(I6,F20.7,ES22.7)") ii,E + E_Nuc,dabs(E - ET)
       end if
       
       if(dabs(Et - E) < 1.0d-9) exit
       
       ET = E
       
    end do
    
    write(*,*) "*****************************************************"
    write(100,*) "*****************************************************"
    
    if( ii <= 9999 ) then
        write(*,"(A17,I6,A9)") " CI Done , After",ii," cycles."
        write(100,"(A17,I6,A9)") " CI Done , After",ii," cycles."
        write(*,*) "*****************************************************"
        write(100,*) "*****************************************************"
    end if 
    
end subroutine Dig_CI_Ham
!=======================================================================
subroutine Read_Mem(n_mem,n_numb)
    ! n_mem : MemFree , n_numb : The max number of floating number.
    implicit none
    integer n_mem,n_numb
    character(len=30) c_mem
   
    open(unit=64,file="/proc/meminfo")
    read(64,*)
    read(64,*) c_mem,n_mem
    close(64)
   
    n_numb = int(real(n_mem)*100d0)
   
end subroutine Read_Mem
!=======================================================================
