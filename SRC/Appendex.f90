subroutine Filename(fname,sfname,rec_char)
    ! Read the file's name,and remove the path.
    implicit none
    integer N,i,j,rec_char
    character(len=200) fname,sfname
    
    N = len(fname)
    j = 0
    rec_char = 0
    
    do i = N,1,-1
         if( fname(i:i) == "/" ) then
             j = i
             exit
         end if
    end do 

    do i = j+1,N
         if( fname(i:i) == "." ) exit
         rec_char = rec_char + 1
         sfname(rec_char:rec_char) = fname(i:i)
    end do
       
end subroutine Filename
!=======================================================================
subroutine Head_Prient(t1,d1,i)
    
    implicit none
    integer i
    character(len=10) t1
    character(len=8) d1
    
    call date_and_time( date = d1 , time = t1)
    
    if( i == 1) then
    write(*,*)
    write(*,*) " ******************************************************** "
    write(*,*) " *                        Amesp                         * "
    write(*,*) " *  Amateurish molecular electronic structure program.  * "
    write(*,*) " *          Warm_Cloud : <3355196386@qq.com>            * "
    write(*,*) " ******************************************************** "
    write(*,*)
    write(*,*) 
    end if
    
    if( i == 2 ) then
    write(100,*) " *************************************************************** "
    write(100,*) " *                           Amesp                             * "
    write(100,*) " *     Amateurish molecular electronic structure program.      * "
    write(100,*) " *             Warm_Cloud : <3355196386@qq.com>                * "
    write(100,*) " *************************************************************** "
    write(100,*)
    write(100,*)
    end if
    
end subroutine Head_Prient
!=======================================================================
subroutine Tail_Prient(t1,d1)

    implicit none
    character(len=10) t1,t2
    character(len=8) d1,d2
    
    call date_and_time( date = d2 , time = t2)
    write(*,*)
    write(*,*)
    write(*,*) "The program run from  ",d1(1:4),"-",d1(5:6),"-", &
                d1(7:8),"  ",t1(1:2),":",t1(3:4),":",t1(5:6),"  to  ",&
                d2(1:4),"-",d2(5:6),"-",d2(7:8),"  ",t2(1:2),":",&
                t2(3:4),":",t2(5:6)
                
    write(100,*)
    write(100,*) "===================================================================="
    write(100,*)
    write(100,*) "The program run from  ",d1(1:4),"-",d1(5:6),"-", &
                d1(7:8),"  ",t1(1:2),":",t1(3:4),":",t1(5:6),"  to  ",&
                d2(1:4),"-",d2(5:6),"-",d2(7:8),"  ",t2(1:2),":",&
                t2(3:4),":",t2(5:6)
                
end subroutine Tail_Prient
!=======================================================================
subroutine Time_Prient(Tim)

    implicit none  
    real(kind=8) tim(4)
      
    write(*,10) " Job cpu time:  ",int(tim(4)),"days,",int(tim(3)),&
                "hours,",int(tim(2)),"minutes,",tim(1),"seconds."
    write(100,10) " Job cpu time:  ",int(tim(4)),"days,",int(tim(3)),&
                "hours,",int(tim(2)),"minutes,",tim(1),"seconds."
                         
    10 format(A15,I5,A6,I4,A8,I4,A9,F8.4,A9) 

end subroutine Time_Prient
!=======================================================================
subroutine Cotime(time,tim)

    implicit none
    integer itime
    real(kind=8) tim(4),time
     
    tim = 0
    itime = int(time)
   
    select case(itime)
   
       case(0:59)
          tim(1) = time   
       case(60:3599)
          tim(2) = itime/60
          tim(1) = time - tim(2)*60
       case(3600:86399)
          tim(3) = itime/3600 
          tim(2) = real(int((time - tim(3)*3600)/60))
          tim(1) = time - tim(3)*3600 - tim(2)*60
       case(86400:2000000000)
          tim(4) = itime/86400
          tim(3) = real(int((time - tim(4)*86400)/3600))
          tim(2) = real(int((time - tim(4)*86400 - tim(3)*3600)/60))
          tim(1) = time - tim(4)*86400 - tim(3)*3600 - tim(2)*60
          
    end select
    
end subroutine
!=======================================================================
subroutine Mat_Out1(A,N0)
    ! Print N*N matrix.
    implicit none  
    integer N0,i,j
    real(kind=8) A(N0,N0)
    character(len=5) Ctemp
    character(len=6) Stn(N0/5*5+1:N0)
      
    do i = 1,N0/5
         write(100,"(I14,4(I13))") 5*i-4,5*i-3,5*i-2,5*i-1,5*i
         do j = 1,N0
              write(100,"(I4,5F13.6)") j,A(j,5*i-4:5*i)
         end do
    end do
    
    if(N0/5*5 == N0) return 
    
    do i = N0/5*5+1,N0
         write(CTemp,"(I4)") i
         Stn(i) = CTemp
    end do
    
    write(100,"(A16,4A13)") Stn(N0/5*5+1:N0)
    
    do j = 1,N0
         write(100,"(I4,5F13.6)") j,A(j,N0/5*5+1:N0)
    end do

end subroutine Mat_Out1
!=======================================================================
subroutine Mat_Out2(A,B,S,N0,IAtom,SAtom,Ind_Basis)
    !Print molecular orbital matrix.
    implicit none  
    integer N0,i,j
    real(kind=8) A(N0,N0),B(N0)
    character(len=1) S(N0)
    character(len=5) Ctemp
    character(len=6) Stn(N0/5*5+1:N0)
    character(len=2) SAtom(N0)
    character(len=4) Ind_Basis(N0)
    character(len=4) IAtom(N0)
       
    do i = 1,N0/5
         write(100,"(I24,4(I11))") 5*i-4,5*i-3,5*i-2,5*i-1,5*i
         write(100,"(A24,4(A11))") S(5*i-4),S(5*i-3),S(5*i-2),S(5*i-1),S(5*i)
         write(100,"(A14,F13.5,1000F11.5)") "Eigenvalues:", B(5*i-4:5*i)
         do j = 1,N0
               write(100,"(I4,A4,A8,5F11.5)") j,IAtom(j),SAtom(j)//" "//Ind_Basis(j),A(j,5*i-4:5*i)
         end do
    end do
    
    if(N0/5*5 == N0) return
    
    do i = N0/5*5+1,N0
         write(CTemp,"(I4)") i
         Stn(i) = CTemp
    end do
    
    write(100,"(A26,4A11)") Stn(N0/5*5+1:N0)
    write(100,"(A24,4(A11))") S(N0/5*5+1:N0)
    write(100,"(A14,F13.5,1000F11.5)") "Eigenvalues:", B(N0/5*5+1:N0)
    
    do j = 1,N0
         write(100,"(I4,A4,A8,5F11.5)") j,IAtom(j),SAtom(j)//" "//Ind_Basis(j),A(j,N0/5*5+1:N0)
    end do
        

end subroutine Mat_Out2
!=======================================================================
subroutine GOP_Out1(Mulliken,N,IAtom,SAtom,Ind_Basis)

    implicit none
    integer N,i
    real(kind=8) Mulliken(N,N)
    character(len=2) SAtom(N) 
    character(len=4) Ind_Basis(N)
    character(len=4) IAtom(N)   

    do i = 1,N
         write(100,"(I4,A4,A8,F14.6)") i,IAtom(i),SAtom(i)//" "//Ind_Basis(i),sum(Mulliken(:,i))
    end do  
    
end subroutine GOP_Out1
!=======================================================================
subroutine GOP_Out2(MullikenA,MullikenB,N,IAtom,SAtom,Ind_Basis)

    implicit none
    integer N,i
    real(kind=8) MullikenA(N,N),MullikenB(N,N)
    character(len=2) SAtom(N) 
    character(len=4) Ind_Basis(N)
    character(len=4) IAtom(N)
    
    write(100,"(A28,3A13)") "Total","Alpha","Beta ","Spin "
    
    do i = 1,N
         write(100,"(I4,A4,A8,4F13.6)") i,IAtom(i),SAtom(i)//" "//Ind_Basis(i),sum(MullikenA(:,i)+MullikenB(:,i)), &
               sum(MullikenA(:,i)),sum(MullikenB(:,i)),sum(MullikenA(:,i)-MullikenB(:,i))
    end do  
    
end subroutine GOP_Out2
!=======================================================================
subroutine Two_Ele_Int_Out(N,Neri,Two_Electron,EI)
    !Print Two_Electron_Int
    implicit none
    integer i,j,k,l,N,Neri,EI(N,N)
    real(kind=8) Two_Electron(N*(N+1)/2,Neri*(Neri+1)/2)

    write(100,*) " Two-electron integral:"
    write(100,*)
    write(100,"(4A6)") "i","j","k","l" 
    write(100,*) 
    do i = 1,N
      do j = i,N
        do k = 1,N
          do l = k,N
             if( i*(i-1)+2*j < k*(k-1)+2*l ) cycle
             if( dabs(Two_Electron(EI(l,k),EI(j,i))) < 1.0d-8 ) cycle
             write(100,"(4I6,F16.8)") i,j,k,l,Two_Electron(EI(l,k),EI(j,i))
          end do
        end do
      end do
    end do
    write(100,*)
          
end subroutine Two_Ele_Int_Out
!=======================================================================
