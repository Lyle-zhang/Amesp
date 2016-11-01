subroutine TDHF_Close(nstate,SCov,En,Ee,g_Mo,EI,MoCu,Dipole,N,Nocc,E_Nuc)

    implicit none 
    integer ii,i,N,Nocc,ntd,EI(N,N),nstate,SCov
    integer icis,jcis,j,a,b
    integer imo,amo,kcis
    integer nwo,nwv,no1,no2,nv1,nv2
    integer,allocatable :: str1(:),str2(:)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2),En(N),E_Nuc,Ee,norm
    real(kind=8) Ptrans(N,N),sqrDtrans,MoCu(N,N),Dipole(N,N,3)
    real(kind=8),allocatable :: A_D(:,:),B_D(:,:),H_D(:,:),VX(:,:),VY(:,:),Dtrans(:,:),osc_str(:)
    real(kind=8),allocatable :: MTemp(:,:),MTemp1(:,:),T_D(:,:),E_D(:),XpY(:,:),XmY(:,:)
    character(len=4) CTemp1,CTemp2
    character(len=13),allocatable :: CX(:),CY(:)
    
    no1 = 1 ; no2 = Nocc ; nv1 = Nocc+1 ; nv2 = N
    nwo = no2 - no1 + 1
    nwv = nv2 - nv1 + 1    
    
    if ( nwv < 1 ) then
        write(*,*)
        write(*,*) " Segmentation Fault!"
        stop
    end if
        
    ntd = nwo*nwv
    
    write(*,*)
    write(*,*) " Ntdhf =", ntd

    write(100,*)"************************************************************"
    write(100,*)"Excited states from <AA,BB:AA,BB> singles matrix:           "
    write(100,*)"************************************************************"  
    write(100,*)
      
    allocate(A_D(ntd,ntd),B_D(ntd,ntd),H_D(ntd,ntd),VX(ntd,ntd),VY(ntd,ntd),str1(ntd), &
             str2(ntd),CX(ntd),CY(ntd),MTemp(ntd,ntd),MTemp1(ntd,ntd),T_D(ntd,ntd), &
             XpY(ntd,ntd),XmY(ntd,ntd),E_D(ntd),Dtrans(3,nstate),osc_str(nstate))
    !----------------------Construct A_D,B_D Matrix---------------------
    kcis = 0
    do imo = no1,no2
       do amo = nv1,nv2
          kcis = kcis + 1
          str1(kcis) = imo
          str2(kcis) = amo
          write(CTemp1,"(i4)") imo
          write(CTemp2,"(i4)") amo
          CX(kcis) = CTemp1//" --> "//CTemp2
          CY(kcis) = CTemp1//" <-- "//CTemp2
       end do
    end do
    
    do icis = 1,ntd
        i = str1(icis)
        a = str2(icis)
        do jcis = icis,ntd       
            j = str1(jcis)
            b = str2(jcis)
            
            A_D(icis,jcis) = 2*g_Mo(EI(i,a),EI(j,b)) - g_Mo(EI(i,j),EI(a,b))
            B_D(icis,jcis) = 2*g_Mo(EI(i,a),EI(j,b)) - g_Mo(EI(i,b),EI(j,a))
            if( i==j .and. a==b ) A_D(icis,jcis) = A_D(icis,jcis) + En(a) - En(i)
            
            A_D(jcis,icis) = A_D(icis,jcis) ! A
            B_D(jcis,icis) = B_D(icis,jcis) ! B
            
        end do
    end do
    !===================================================================
    MTemp = A_D - B_D                         ! (A-B)
    call  Mat_Sqrt1(MTemp,MTemp1,ntd)         ! MTemp1 = (A-B)^(1/2) 
    MTemp = A_D + B_D                         ! (A+B)
    H_D = matmul(matmul(MTemp1,MTemp),MTemp1) ! H = (A-B)^(1/2)*(A+B)*(A-B)^(1/2)
    call SY_Eigen_System(H_D,E_D,T_D,ntd)     ! HT = w^2T
    
    XpY = matmul(MTemp1,T_D)                  ! |X+Y> = (A-B)^(1/2)*T
    E_D = dsqrt(dabs(E_D))                    ! E = w^2,w = sqrt[E])
    
    do i = 1,ntd
         XmY(:,i) = matmul(MTemp,XpY(:,i))/E_D(i) ! w*|X-Y> = (A+B)*|X+Y>
    end do 
    
    VX=(XpY+XmY)/2.0d0                        !  X = ( |X+Y> + |X-Y> )/2
    VY=(XpY-XmY)/2.0d0                        !  Y = ( |X+Y> - |X-Y> )/2  
    
    do i = 1,ntd
        norm=dsqrt(sum(VX(:,i)**2) - sum(VY(:,i)**2))
        !Normalization, meantime convert CSF coefficient to Slater determinant coefficient
        VX(:,i) = VX(:,i)/norm/dsqrt(2.0d0)
        VY(:,i) = VY(:,i)/norm/dsqrt(2.0d0)
    end do 
    !===================== oscillator strength =========================
    do ii = 1,nstate
        Ptrans = 0.0d0
        do icis = 1,ntd
            i = str1(icis)
            j = str2(icis)
            Ptrans =  Ptrans + VX(icis,ii)*matmul(MoCu(:,i:i),transpose(MoCu(:,j:j)))
            Ptrans =  Ptrans + VY(icis,ii)*matmul(MoCu(:,i:i),transpose(MoCu(:,j:j)))
        end do
        Ptrans = Ptrans*2.0d0 !Since we need consider both Alpha->Alpha and Beta->Beta determinants
        Dtrans(1,ii) = sum(Ptrans*Dipole(:,:,1))
        Dtrans(2,ii) = sum(Ptrans*Dipole(:,:,2))
        Dtrans(3,ii) = sum(Ptrans*Dipole(:,:,3))
        sqrDtrans=sum(Dtrans(:,ii)**2)
        osc_str(ii)=2.0d0/3.0d0*E_D(ii)*sqrDtrans
    end do
    !===================================================================
    do ii = 1,nstate    
        write(100,"(A14,I4,A2,F10.4,A3,A2,F10.2,A3,A2,A4,F8.4)") " Excited State",ii," ,",E_D(ii)*27.2114d0, &
                  " ev"," ,",45.5633d0/E_D(ii)," nm"," ,"," f= ",osc_str(ii)
                                                    
        do i = 1,ntd
            if( dabs(VX(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") CX(i),VX(i,ii)
        end do
        
        do i = 1,ntd
            if( dabs(VY(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") CY(i),VY(i,ii)
        end do
        
        write(100,"(A12,F15.8,A3,A10,F7.4)") "E(TD-HF) = ",E_D(ii) + Ee + E_Nuc,"  ,"," <S**2> =",0.0d0
        write(100,*)
    end do
          
    write(100,*)"************************************************************"    

end subroutine TDHF_Close
!=======================================================================
!=======================================================================
subroutine TDHF_Open(nstate,SCov,EnA,EnB,Ee,g_MoA,g_MoB,g_MoAB,EI,MoCuA, &
                     MoCuB,Dipole,N,NoccA,NoccB,E_Nuc)

    implicit none 
    integer nstate,SCov,ii,N,NoccA,NoccB,ntd,EI(N,N)
    integer i,j,a,b,imo,amo,kcis,icis,jcis
    integer nwoA,nwvA,no1A,no2A,nv1A,nv2A
    integer nwoB,nwvB,no1B,no2B,nv1B,nv2B
    integer,allocatable :: str1(:),str2(:)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) PtransA(N,N),PtransB(N,N),sqrDtrans,MoCuA(N,N),MoCuB(N,N),Dipole(N,N,3)
    real(kind=8) EnA(N),EnB(N),Ee,E_Nuc,norm,Ptrans(N,N)
    real(kind=8),allocatable :: A_D(:,:),B_D(:,:),H_D(:,:),VX(:,:),VY(:,:),Dtrans(:,:),osc_str(:)
    real(kind=8),allocatable :: MTemp(:,:),MTemp1(:,:),T_D(:,:),E_D(:),XpY(:,:),XmY(:,:)
    character(len=13),allocatable :: CX(:),CY(:)
    character(len=4) CTemp1,CTemp2
    
    no1A = 1 ; no2A = NoccA ; nv1A = NoccA+1 ; nv2A = N
    nwoA = no2A - no1A + 1
    nwvA = nv2A - nv1A + 1
    
    no1B = 1 ; no2B = NoccB ; nv1B = NoccB+1 ; nv2B = N
    nwoB = no2B - no1B + 1
    nwvB = nv2B - nv1B + 1    

    if ( nwvA < 1 .or. nwvB < 1 ) then
        write(*,*)
        write(*,*) " Segmentation Fault!"
        stop
    end if
       
    ntd = nwoA*nwvA + nwoB*nwvB
    write(*,*)
    write(*,*) " Ntdhf =", ntd
    
    write(100,*)"************************************************************"
    write(100,*)"Excited states from <AA,BB:AA,BB> singles matrix:           "
    write(100,*)"************************************************************"  
    write(100,*)
    
    allocate(A_D(ntd,ntd),B_D(ntd,ntd),H_D(ntd,ntd),VX(ntd,ntd),VY(ntd,ntd),str1(ntd), &
             str2(ntd),CX(ntd),CY(ntd),MTemp(ntd,ntd),MTemp1(ntd,ntd),T_D(ntd,ntd), &
             XpY(ntd,ntd),XmY(ntd,ntd),E_D(ntd),Dtrans(3,nstate),osc_str(nstate))
             
    !=====================Construct A_D,B_D Matrix-=====================
    kcis = 0
    do imo = no1A,no2A
       do amo = nv1A,nv2A
          kcis = kcis + 1
          str1(kcis) = imo
          str2(kcis) = amo
          write(CTemp1,"(i4)") imo
          write(CTemp2,"(i4)") amo
          CX(kcis) = CTemp1//" --> "//CTemp2
          CY(kcis) = CTemp1//" <-- "//CTemp2
       end do
    end do
    
    do imo = no1B,no2B
       do amo = nv1B,nv2B
          kcis = kcis + 1
          str1(kcis) = imo
          str2(kcis) = amo
          write(CTemp1,"(i4)") imo
          write(CTemp2,"(i4)") amo
          CX(kcis) = CTemp1//" --> "//CTemp2
          CY(kcis) = CTemp1//" <-- "//CTemp2
       end do
    end do
    !---------------------------alpha-alpha-----------------------------
    do icis = 1,nwoA*nwvA
        i = str1(icis)
        a = str2(icis)
        do jcis = icis,nwoA*nwvA      
            j = str1(jcis)
            b = str2(jcis)
            
            A_D(icis,jcis) = g_MoA(EI(i,a),EI(j,b)) - g_MoA(EI(i,j),EI(a,b))
            B_D(icis,jcis) = g_MoA(EI(i,a),EI(j,b)) - g_MoA(EI(i,b),EI(j,a))
            if( i==j .and. a==b ) A_D(icis,jcis) = A_D(icis,jcis) + EnA(a) - EnA(i)
            
        end do
    end do
    !----------------------------beta-beta------------------------------
    do icis = nwoA*nwvA+1,ntd
        i = str1(icis)
        a = str2(icis)
        do jcis = icis,ntd     
            j = str1(jcis)
            b = str2(jcis)
            
            A_D(icis,jcis) = g_MoB(EI(i,a),EI(j,b)) - g_MoB(EI(i,j),EI(a,b))
            B_D(icis,jcis) = g_MoB(EI(i,a),EI(j,b)) - g_MoB(EI(i,b),EI(j,a))
            if( i==j .and. a==b ) A_D(icis,jcis) = A_D(icis,jcis) + EnB(a) - EnB(i)
            
        end do
    end do
    !----------------------------alpha-beta-----------------------------
    do icis = 1,nwoA*nwvA
        i = str1(icis)
        a = str2(icis)
        do jcis = nwoA*nwvA+1,ntd
            j = str1(jcis)
            b = str2(jcis)
            
            A_D(icis,jcis) = g_MoAB(EI(j,b),EI(i,a))
            B_D(icis,jcis) = g_MoAB(EI(j,b),EI(i,a))
            
        end do
    end do
    !-------------------------------------------------------------------
    do icis = 1,ntd
        do jcis = icis,ntd
            A_D(jcis,icis) = A_D(icis,jcis)
            B_D(jcis,icis) = B_D(icis,jcis)
        end do
    end do
    !===================================================================
    MTemp = A_D - B_D                         ! (A-B)
    call  Mat_Sqrt1(MTemp,MTemp1,ntd)         ! MTemp1 = (A-B)^(1/2) 
    MTemp = A_D + B_D                         ! (A+B)
    H_D = matmul(matmul(MTemp1,MTemp),MTemp1) ! H = (A-B)^(1/2)*(A+B)*(A-B)^(1/2)
    call SY_Eigen_System(H_D,E_D,T_D,ntd)     ! HT = w^2T
    
    XpY = matmul(MTemp1,T_D)                  ! |X+Y> = (A-B)^(1/2)*T
    E_D = dsqrt(dabs(E_D))                    ! E = w^2,w = sqrt[E])
    
    do i = 1,ntd
         XmY(:,i) = matmul(MTemp,XpY(:,i))/E_D(i) ! w*|X-Y> = (A+B)*|X+Y>
    end do 
    
    VX=(XpY+XmY)/2.0d0                        !  X = ( |X+Y> + |X-Y> )/2
    VY=(XpY-XmY)/2.0d0                        !  Y = ( |X+Y> - |X-Y> )/2  

    do i = 1,ntd
        norm=dsqrt(sum(VX(:,i)**2) - sum(VY(:,i)**2))
        VX(:,i) = VX(:,i)/norm
        VY(:,i) = VY(:,i)/norm
    end do 
    !===================== oscillator strength =========================
    do ii = 1,nstate
        PtransA = 0.0d0
        PtransB = 0.0d0
        
        do icis = 1,nwoA*nwvA
            i = str1(icis)
            j = str2(icis)
            PtransA =  PtransA + VX(icis,ii)*matmul(MoCuA(:,i:i),transpose(MoCuA(:,j:j)))
            PtransA =  PtransA + VY(icis,ii)*matmul(MoCuA(:,i:i),transpose(MoCuA(:,j:j)))
        end do
        
        do icis = nwoA*nwvA + 1,ntd
            i = str1(icis)
            j = str2(icis)
            PtransB =  PtransB + VX(icis,ii)*matmul(MoCuB(:,i:i),transpose(MoCuB(:,j:j)))
            PtransB =  PtransB + VY(icis,ii)*matmul(MoCuB(:,i:i),transpose(MoCuB(:,j:j)))
        end do
        
        Ptrans = PtransA + PtransB
        Dtrans(1,ii) = sum(Ptrans*Dipole(:,:,1))
        Dtrans(2,ii) = sum(Ptrans*Dipole(:,:,2))
        Dtrans(3,ii) = sum(Ptrans*Dipole(:,:,3))
        sqrDtrans=sum(Dtrans(:,ii)**2)
        osc_str(ii)=2.0d0/3.0d0*E_D(ii)*sqrDtrans
    end do
    !===================================================================
    
    do ii = 1,nstate
    
        write(100,"(A14,I4,A2,F10.4,A3,A2,F10.2,A3,A2,A4,F8.4)") " Excited State",ii," ,",E_D(ii)*27.2114d0, &
                    " ev"," ,",45.5633d0/E_D(ii)," nm"," ,"," f= ",osc_str(ii)
        write(100,*) "    Alpha State:"
        do i = 1,nwoA*nwvA
             if( dabs(VX(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") CX(i),VX(i,ii)
        end do
        do i = 1,nwoA*nwvA
             if( dabs(VY(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") CY(i),VY(i,ii)
        end do
            
        write(100,*) "    Beta State:"
        do i = nwoA*nwvA + 1,ntd
             if( dabs(VX(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") CX(i),VX(i,ii)
        end do
        do i = nwoA*nwvA + 1,ntd
             if( dabs(VY(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") CY(i),VY(i,ii)
        end do
        
        write(100,"(A12,F15.8)") "E(TD-HF) = ",E_D(ii) + Ee + E_Nuc
        write(100,*)  
        
    end do
         
    write(100,*)"************************************************************"      
      

end subroutine TDHF_Open
!=======================================================================


