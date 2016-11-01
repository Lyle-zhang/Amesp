subroutine Electron_Density(N,Density,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
    
    implicit none
    integer,parameter :: N0=100
    integer i,j,k,l,m,N,l1,ii,jj,bas_num(N),NAtom(N),l_Basis(3,N)
    real(kind=8) x1,x2,y1,y2,z1,z2,hx,hy,hz
    real(kind=8) xt1(N0*N0*N0),yt1(N0*N0*N0),zt1(N0*N0*N0)
    real(kind=8) x(N0),y(N0),z(N0),Xi(N),rho(N0,N0,N0),Temp1,BN
    real(kind=8) Density(N,N),c_Basis(6,N,36),a_Basis(6,N,36),Rt(3,N)  
    character(len=1) ch
    logical log1
    
    x2 = maxval(Rt(1,:)) + 2.0d0
    x1 = minval(Rt(1,:)) - 2.0d0
    y2 = maxval(Rt(2,:)) + 2.0d0
    y1 = minval(Rt(2,:)) - 2.0d0
    z2 = maxval(Rt(3,:)) + 2.0d0
    z1 = minval(Rt(3,:)) - 2.0d0
    hx = (x2-x1)/real(N0,8)
    hy = (y2-y1)/real(N0,8)
    hz = (z2-z1)/real(N0,8)

    write(*,*)
    write(*,*) " Press y to draw (spin) density map,or c to continue:"
    read(*,*) ch
    if(ch /= "y") return
              
    do i = 1,N0
      do j = 1,N0
        do m = 1,N0
         
           x(i) = x1 + i*hx
           y(j) = y1 + j*hy
           z(m) = z1 + m*hz
           
           do ii = 1,N
         
              Xi(ii) = 0.0d0
              do jj = 1,bas_num(ii)
              call Bas_Normal(a_Basis(jj,ii,NAtom(ii)),l_Basis(1,ii),l_Basis(2,ii),l_Basis(3,ii),BN)
              Xi(ii) = c_Basis(jj,ii,NAtom(ii))*BN* &
                       dexp(-a_Basis(jj,ii,NAtom(ii))*((x(i)-Rt(1,ii))**2+ &
                       (y(j)-Rt(2,ii))**2+(z(m)-Rt(3,ii))**2)) + Xi(ii)
              end do
              Xi(ii) = Xi(ii)*(x(i)-Rt(1,ii))**l_Basis(1,ii)*(y(j)-Rt(2,ii))**l_Basis(2,ii) &
                       *(z(m)-Rt(3,ii))**l_Basis(3,ii)                  
           end do 
         
           rho(m,j,i) = 0.0d0         
           do k = 1,N
              Temp1 = 0.0d0
              do l = 1,N
                 Temp1 = Temp1 + Density(k,l)*Xi(k)*Xi(l)
              end do
              rho(m,j,i) = rho(m,j,i) + Temp1
           end do

        end do
      end do
    end do
    
    if (ch == "y") then
       l1 = 0 
       do i = 1,N0
         do j = 1,N0
           do m = 1,N0
              log1 = rho(m,j,i) > 0.08d0 .and. rho(m,j,i) < 0.12d0 
              if(log1) then
                  l1 = l1 + 1
                  xt1(l1) = x(i)
                  yt1(l1) = y(j)
                  zt1(l1) = z(m)
              end if
            end do
         end do
       end do
       call Plot3D(xt1,yt1,zt1,l1)
       call System("rm Data.txt")
    end if
    
             
end subroutine Electron_Density
!=======================================================================
subroutine Molecular_Orbital(N,MoCu,En,c_Basis,a_Basis,l_Basis,Rt,NAtom,bas_num)
    
    implicit none
    integer,parameter :: N0=100
    integer i,j,k,l,m,N,ii,jj,bas_num(N),NAtom(N),l1,l2,l_Basis(3,N)
    real(kind=8) x1,x2,y1,y2,z1,z2,hx,hy,hz,x(N0),y(N0),z(N0),BN
    real(kind=8) xt1(N0*N0*N0),yt1(N0*N0*N0),zt1(N0*N0*N0),Phi(N0,N0,N0),Xi(N),Temp1
    real(kind=8) xt2(N0*N0*N0),yt2(N0*N0*N0),zt2(N0*N0*N0)
    real(kind=8) MoCu(N,N),En(N),c_Basis(6,N,36),a_Basis(6,N,36),Rt(3,N)  
    character(len=1) ch
    logical log1,log2
    

    x2 = maxval(Rt(1,:)) + 2.0d0
    x1 = minval(Rt(1,:)) - 2.0d0
    y2 = maxval(Rt(2,:)) + 2.0d0
    y1 = minval(Rt(2,:)) - 2.0d0
    z2 = maxval(Rt(3,:)) + 2.0d0
    z1 = minval(Rt(3,:)) - 2.0d0
    hx = (x2-x1)/real(N0,8)
    hy = (y2-y1)/real(N0,8)
    hz = (z2-z1)/real(N0,8)
    
101 write(*,*)
    write(*,*) " Which orbital you want to draw ?"
    write(*,"(A15,I5)") " input k , k<=",N 
    read(*,*) k
    if( k > N ) stop
    write(*,"(A4,G0,A3,ES20.8)") " E(",k,") =",En(k)   
    !-------------------------------------------------------------------   
    do i = 1,N0
      do j = 1,N0
         do m = 1,N0
         
           x(i) = x1 + i*hx
           y(j) = y1 + j*hy
           z(m) = z1 + m*hz
           do ii = 1,N
         
              Xi(ii) = 0.0d0
              do jj = 1,bas_num(ii)
              call Bas_Normal(a_Basis(jj,ii,NAtom(ii)),l_Basis(1,ii),l_Basis(2,ii),l_Basis(3,ii),BN)
              Xi(ii) = c_Basis(jj,ii,NAtom(ii))*BN* &
                       dexp(-a_Basis(jj,ii,NAtom(ii))*((x(i)-Rt(1,ii))**2+ &
                       (y(j)-Rt(2,ii))**2+(z(m)-Rt(3,ii))**2)) + Xi(ii)
              end do
              Xi(ii) = Xi(ii)*(x(i)-Rt(1,ii))**l_Basis(1,ii)*(y(j)-Rt(2,ii))**l_Basis(2,ii) &
                       *(z(m)-Rt(3,ii))**l_Basis(3,ii)                  
           end do 
         
           Temp1 = 0.0d0
           do l = 1,N
               Temp1 = Temp1 + MoCu(l,k)*Xi(l)
           end do

           phi(m,j,i) = Temp1
         
         end do
      end do
    end do
    !-------------------------------------------------------------------
    l1 = 0 ; l2 = 0
    do i = 1,N0
      do j = 1,N0
         do m = 1,N0
             log1 = phi(m,j,i) > 4.0d-2 .and. phi(m,j,i) < 4.2d-2 
             log2 = phi(m,j,i) < -4.0d-2 .and. phi(m,j,i) > -4.2d-2 
             if(log1) then
                 l1 = l1 + 1
                 xt1(l1) = x(i)
                 yt1(l1) = y(j)
                 zt1(l1) = z(m)
             end if
             if(log2) then
                 l2 = l2 + 1
                 xt2(l2) = x(i)
                 yt2(l2) = y(j)
                 zt2(l2) = z(m)
             end if
         end do
      end do
    end do
    !-------------------------------------------------------------------
    call surf(xt1,yt1,zt1,l1,xt2,yt2,zt2,l2)
    call System("rm Data1.txt Data2.txt")
    
    write(*,*) " Press c to contine,or q to quit:"
    read(*,*) ch
    
    if (ch == "c") goto 101
    
             
end subroutine Molecular_orbital
!=======================================================================
subroutine Plot3D(x,y,z,N0)

    implicit none
    integer i,N0
    integer status
    integer system
    real(kind=8) x(N0),y(N0),z(N0)
    character(len = 200) Command
    
    open( unit=300,file="Data.txt" )
    open( unit=301,file="/tmp/Command.txt" )
    
    do i = 1,N0
       write(300,"(3ES20.9)") x(i),y(i),z(i)
    end do
    write(300,*)
    close(300)
       
    write(301,*) 'set terminal wxt  title  "Gnuplot"'
    write(301,*) "set nokey"
    write(301,*) 'splot "Data.txt" '
    write(301,*) 'pause -1 "  Press Enter to continue ..."'
    write(301,*) "q"
    close(301)

    write (Command, *) 'gnuplot '//trim("/tmp/Command.txt")
    status = system(trim(Command))
	
    if (status /= 0) then
        write(*,*) " Fail to run gnuplot!"
        stop
    end if
	
end subroutine Plot3D
!=======================================================================
subroutine Surf(x1,y1,z1,N1,x2,y2,z2,N2)

    implicit none
    integer i,N1,N2
    integer status
    integer system
    real(kind=8) x1(N1),y1(N1),z1(N1)
    real(kind=8) x2(N2),y2(N2),z2(N2)
    character(len = 200) Command
    
    open( unit=300,file="Data1.txt" )
    open( unit=302,file="Data2.txt" )
    open( unit=301,file="/tmp/Command.txt" )
    
    do i = 1,N1
       write(300,"(3ES20.9)") x1(i),y1(i),z1(i)
    end do
    write(300,*)
    close(300)
    
    do i = 1,N2
       write(302,"(3ES20.9)") x2(i),y2(i),z2(i)
    end do
    write(302,*)
    close(302)
       
    write(301,*) 'set terminal wxt  title  "Gnuplot"'
    write(301,*) "set nokey"
    write(301,"(A82)") ' splot "Data1.txt" with points lc rgb "red","Data2.txt" with points lc rgb "green"'
    write(301,*) 'pause -1 "  Press Enter to continue ..."'
    write(301,*) "q"
    close(301)

    write (Command, *) 'gnuplot '//trim("/tmp/Command.txt")
    status = system(trim(Command))
	
    if (status /= 0) then
        write(*,*) " Fail to run gnuplot!"
        stop
    end if
	
end subroutine Surf
!=======================================================================
