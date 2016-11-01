subroutine CIS_Close(nstate,SCov,En,Ee,g_Mo,EI,MoCu,Dipole,N,Nocc,E_Nuc)

    implicit none 
    integer ii,i,N,Nocc,ncis,EI(N,N),nstate,SCov
    integer icis,jcis,j,a,b
    integer imo,amo,kcis
    integer nwo,nwv,no1,no2,nv1,nv2
    integer,allocatable :: str1(:),str2(:)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2),En(N),E_Nuc,Ee
    real(kind=8) Ptrans(N,N),sqrDtrans,MoCu(N,N),Dipole(N,N,3)
    real(kind=8),allocatable :: CIS_Ham(:,:),Ecis(:),Vcis(:,:),Dtrans(:,:),osc_str(:)
    character(len=4) CTemp1,CTemp2
    character(len=13),allocatable :: Ccis(:)
    
    no1 = 1 ; no2 = Nocc ; nv1 = Nocc+1 ; nv2 = N
    nwo = no2 - no1 + 1
    nwv = nv2 - nv1 + 1    
    
    if ( nwv < 1 ) then
        write(*,*)
        write(*,*) " Segmentation Fault!"
        stop
    end if
        
    ncis = nwo*nwv
    
    write(*,*)
    write(*,*) " Ncis =", ncis 

    write(100,*)"************************************************************"
    write(100,*)"Excited states from <AA,BB:AA,BB> singles matrix:           "
    write(100,*)"************************************************************"  
    write(100,*)
      
    allocate(CIS_Ham(ncis,ncis),Ecis(ncis),Vcis(ncis,ncis),str1(ncis), &
             str2(ncis),Ccis(ncis),Dtrans(3,nstate),osc_str(nstate))
    !-----------------------Construct CIS Matrix------------------------
    kcis = 0
    do imo = no1,no2
       do amo = nv1,nv2
          kcis = kcis + 1
          str1(kcis) = imo
          str2(kcis) = amo
          write(CTemp1,"(i4)") imo
          write(CTemp2,"(i4)") amo
          Ccis(kcis) = CTemp1//" --> "//CTemp2
       end do
    end do
    
    do icis = 1,ncis
        i = str1(icis)
        a = str2(icis)
        do jcis = icis,ncis       
            j = str1(jcis)
            b = str2(jcis)
            
            CIS_Ham(icis,jcis) = 2*g_Mo(EI(i,a),EI(j,b)) - g_Mo(EI(i,j),EI(a,b))
            if( i==j .and. a==b ) CIS_Ham(icis,jcis) = CIS_Ham(icis,jcis) + En(a) - En(i)
            
            CIS_Ham(jcis,icis) = CIS_Ham(icis,jcis)
            
        end do
    end do
    !-------------------------------------------------------------------
    call SY_Eigen_System(CIS_Ham,Ecis,Vcis,ncis)
    
    !--------------------- oscillator strength -------------------------
    do ii = 1,nstate
        Ptrans = 0.0d0
        Vcis(:,ii) = Vcis(:,ii)/dsqrt(2.0d0)
        do icis = 1,ncis
            i = str1(icis)
            j = str2(icis)
            Ptrans =  Ptrans + Vcis(icis,ii)*matmul(MoCu(:,i:i),transpose(MoCu(:,j:j)))
        end do
        Ptrans = Ptrans*2.0d0 !Since we need consider both Alpha->Alpha and Beta->Beta determinants
        Dtrans(1,ii) = sum(Ptrans*Dipole(:,:,1))
        Dtrans(2,ii) = sum(Ptrans*Dipole(:,:,2))
        Dtrans(3,ii) = sum(Ptrans*Dipole(:,:,3))
        sqrDtrans=sum(Dtrans(:,ii)**2)
        osc_str(ii)=2.0d0/3.0d0*Ecis(ii)*sqrDtrans
    end do
    !-------------------------------------------------------------------
    
    do ii = 1,nstate    
        write(100,"(A14,I4,A2,F10.4,A3,A2,F10.2,A3,A2,A4,F8.4)") " Excited State",ii," ,",Ecis(ii)*27.2114d0, &
                               " ev"," ,",45.5633d0/Ecis(ii)," nm"," ,"," f= ",osc_str(ii)
        do i = 1,nwo*nwv
            if( dabs(Vcis(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") Ccis(i),Vcis(i,ii)
        end do
        write(100,"(A10,F15.8,A3,A10,F7.4)") "E(cis) = ",Ecis(ii) + Ee + E_Nuc,"  ,"," <S**2> =",0.0d0
        write(100,*)
    end do
        
    write(100,*)"************************************************************"    

end subroutine CIS_Close
!=======================================================================
subroutine CID_Close(f_Mo,g_Mo,EI,N,Nocc,E_Nuc,Ecid)

    implicit none 
    integer N,Nocc,ncsf,mcsf,EI(N,N)
    integer nwo,nwv,no1,no2,nv1,nv2
    integer n_mem,n_numb
    integer,allocatable :: stra(:,:),strb(:,:)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2),Ecid,f_Mo(N,N),E_Nuc
    real(kind=8),allocatable :: CID_Ham(:,:)
    
    no1 = 1 ; no2 = Nocc ; nv1 = Nocc+1 ; nv2 = N
    nwo = no2 - no1 + 1
    nwv = nv2 - nv1 + 1    
    
    if ( nwv < 1 ) then
        write(*,*)
        write(*,*) " Segmentation Fault!"
        stop
    end if
        
    call Num_Csf_Close(nwo,nwv,ncsf,mcsf)
    call Read_Mem(n_mem,n_numb)
    
    if ( sqrt(real(n_numb)) < real(mcsf) ) then
        write(*,*)
        write(*,*) " Ncid =", mcsf 
        write(*,*)
        write(*,*) " Error : insufficient virtual memory for cid!"
        write(100,*) " Error : insufficient virtual memory for cid!"
        Ecid = - E_Nuc
        return
    end if
    
    write(*,*)
    write(*,*) " Ncid =", mcsf 
    
    allocate(CID_Ham(mcsf,mcsf),stra(no2,ncsf), &
             strb(no2,ncsf))
    call CI_State_Close(stra,strb,ncsf,no1,no2,nv1,nv2)
    call CISD_Hamilton_Close(N,CID_Ham,mcsf,no2,stra,strb,g_Mo,EI,f_Mo)  
      
    call Dig_CI_Ham(CID_Ham,Ecid,mcsf,E_Nuc)

end subroutine CID_Close
!=======================================================================
subroutine CISD_Close(f_Mo,g_Mo,EI,N,Nocc,E_Nuc,Ecisd)

    implicit none 
    integer N,Nocc,ncsf,mcsf,EI(N,N)
    integer nwo,nwv,no1,no2,nv1,nv2
    integer n_mem,n_numb
    integer,allocatable :: stra(:,:),strb(:,:)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2),Ecisd,f_Mo(N,N),E_Nuc
    real(kind=8),allocatable :: CISD_Ham(:,:)
    
    no1 = 1 ; no2 = Nocc ; nv1 = Nocc+1 ; nv2 = N
    nwo = no2 - no1 + 1
    nwv = nv2 - nv1 + 1 

    if ( nwv < 1 ) then
        write(*,*)
        write(*,*) " Segmentation Fault!"
        stop
    end if
        
    call Num_Csf_Close(nwo,nwv,ncsf,mcsf)   
    call Read_Mem(n_mem,n_numb)
    
    if ( sqrt(real(n_numb)) < real(mcsf) ) then
        write(*,*)
        write(*,*) " Ncisd =", ncsf 
        write(*,*)
        write(*,*) " Error : insufficient virtual memory for cisd!"
        write(100,*) " Error : insufficient virtual memory for cisd!"
        Ecisd = - E_Nuc
        return
    end if
    
    write(*,*)
    write(*,*) " Ncisd =", ncsf  
    
    allocate(CISD_Ham(ncsf,ncsf),stra(no2,ncsf), &
             strb(no2,ncsf))
    call CI_State_Close(stra,strb,ncsf,no1,no2,nv1,nv2)
    call CISD_Hamilton_Close(N,CISD_Ham,ncsf,no2,stra,strb,g_Mo,EI,f_Mo) 
         
    call Dig_CI_Ham(CISD_Ham,Ecisd,ncsf,E_Nuc)

end subroutine CISD_Close
!=======================================================================
subroutine CISD_Hamilton_Close(N,CI_Ham,ncsf,no2,stra,strb,g_Mo,EI,f_Mo)

    implicit none
    integer N,i,j,kcsf1,kcsf2,ncsf,no2,ndif,ndifa,ndifb,ipa,ipb
    integer iobt(no2*2),p,q,r,s,EI(N,N)
    integer stra(no2,ncsf),strb(no2,ncsf)
    integer strxa(no2),strya(no2),strma(no2),strna(no2)
    integer strxb(no2),stryb(no2),strmb(no2),strnb(no2)
    real(kind=8) g_Mo(N*(N+1)/2,N*(N+1)/2),f_Mo(N,N)
    real(kind=8) CI_Ham(ncsf,ncsf)
    real(kind=8) Temp1,Temp2,Temp3,Temp4
    
    strxa = 0 ; strya = 0 ; strma = 0 ; strna = 0
    strxb = 0 ; stryb = 0 ; strmb = 0 ; strnb = 0
    CI_Ham = 0.0d0
    
    do kcsf1 = 1,ncsf
       do kcsf2 = kcsf1,ncsf
       
           call Match(no2,stra(:,kcsf1),stra(:,kcsf2),strxa,strya, &
                      strma,strna,ndifa,ipa)
           call Match(no2,strb(:,kcsf1),strb(:,kcsf2),strxb,stryb, &
                      strmb,strnb,ndifb,ipb)
           ndif = ndifa + ndifb

           if ( ndif > 2 ) cycle
           !------------------------------------------------------------
           if ( ndif == 2 ) then

               if ( ndifa == 2 ) then
                   p = strxa(1) ; q = strxa(2)
                   r = strya(1) ; s = strya(2)
                   CI_Ham(kcsf2,kcsf1) = g_Mo(EI(s,q),EI(r,p))-g_Mo(EI(r,q),EI(s,p))
                   
               end if
               
               if ( ndifb == 2 ) then
                   p = strxb(1) ; q = strxb(2)
                   r = stryb(1) ; s = stryb(2)
                   CI_Ham(kcsf2,kcsf1) = g_Mo(EI(s,q),EI(r,p))-g_Mo(EI(r,q),EI(s,p))
                   
               end if
               
               if ( ndifa == 1 ) then
                   p = strxa(1) ; q = strxb(1)
                   r = strya(1) ; s = stryb(1)
                   CI_Ham(kcsf2,kcsf1) = g_Mo(EI(s,q),EI(r,p))
                                     
               end if
               CI_Ham(kcsf2,kcsf1) = CI_Ham(kcsf2,kcsf1)*ipa*ipb
           end if
           !------------------------------------------------------------
           if ( ndif == 1 ) then
               if( kcsf1 == 1 .or. kcsf2 == 1 ) cycle
               
               if( strxa(1) == 0 ) then   
                   iobt(1:no2) = stra(:,kcsf1) 
                   iobt(no2+1:2*no2) = strb(:,kcsf1) 
                   p = strxb(1);q = stryb(1)
                   
                   Temp1 = 0.0d0
                   do j = 1,2*no2 
                                    
                        i = iobt(j) 
                        if( j <= no2 ) then
                            Temp2 = g_Mo(EI(i,i),EI(p,q))
                        else
                            Temp2 = g_Mo(EI(i,i),EI(p,q))-g_Mo(EI(q,i),EI(i,p))
                        end if
                         
                        Temp1 = Temp1 + Temp2
                   end do

                   CI_Ham(kcsf2,kcsf1) = f_Mo(q,p) + Temp1

               end if
               
               if( strxb(1) == 0 ) then   
                   iobt(1:no2) = stra(:,kcsf1) 
                   iobt(no2+1:2*no2) = strb(:,kcsf1) 
                   p = strxa(1);q = strya(1)
                   
                   Temp1 = 0.0d0
                   do j = 1,2*no2  
                                    
                        i = iobt(j) 
                        if( j > no2 ) then
                            Temp2 = g_Mo(EI(i,i),EI(p,q))
                        else
                            Temp2 = g_Mo(EI(i,i),EI(p,q))-g_Mo(EI(q,i),EI(i,p))
                        end if
                         
                        Temp1 = Temp1 + Temp2
                   end do

                   CI_Ham(kcsf2,kcsf1) = f_Mo(q,p) + Temp1

               end if
               
               CI_Ham(kcsf2,kcsf1) = CI_Ham(kcsf2,kcsf1)*ipa*ipb
           end if
           !------------------------------------------------------------
           if ( ndif == 0 ) then
           
                iobt(1:no2) = stra(:,kcsf1) ; iobt(no2+1:2*no2) = strb(:,kcsf1)
               
                Temp1 = 0.0d0
                do i = 1,2*no2
                     p = iobt(i)
                     Temp1 = Temp1 + f_Mo(p,p)
                end do
               
                Temp2 = 0.0d0
                do i = 1,2*no2 
                    p = iobt(i)
                    Temp3 = 0.0d0
                    do j = 1,2*no2                  
                         q = iobt(j) 
                         if( i == j ) cycle 

                         if( j <= no2 .and. i > no2 ) then
                             Temp4 = 0.5d0*g_Mo(EI(q,q),EI(p,p))
                         else if( i <= no2 .and. j > no2 ) then
                             Temp4 = 0.5d0*g_Mo(EI(q,q),EI(p,p)) 
                         else
                             Temp4 = 0.5d0*(g_Mo(EI(q,q),EI(p,p))-g_Mo(EI(q,p),EI(p,q)))
                         end if
                         
                         Temp3 = Temp3 + Temp4
                    end do
                    Temp2 = Temp2 + Temp3
                end do
               
                CI_Ham(kcsf2,kcsf1) = Temp1 + Temp2
           end if   
           CI_Ham(kcsf1,kcsf2) = CI_Ham(kcsf2,kcsf1)    
       end do
    end do
    
    write(*,*) " CI Hamilton complete."
                   
end subroutine CISD_Hamilton_Close
!=======================================================================
subroutine CI_State_Close(stra,strb,ncsf,no1,no2,nv1,nv2)

    implicit none
    integer imo,jmo,amo,bmo,ncsf,no1,no2,nv1,nv2,kcsf
    integer ipa,ipb
    integer stra(no2,ncsf),strb(no2,ncsf)
    
    !--------Ground state--------
    do imo = 1,no2
       stra(imo,:) = imo
       strb(imo,:) = imo
    end do
    
    kcsf = 1 
    
    if(no1 == no2 .or. nv1 == nv2 ) goto 1000
    !--------- a2/b0 .and. a0/b2 ----------
    do imo = no1+1,no2
       do jmo = no1,imo-1
          do amo = nv1+1,nv2
             do bmo = nv1,amo-1
                kcsf = kcsf + 1
                stra(imo,kcsf) = amo
                stra(jmo,kcsf) = bmo
             end do
          end do
       end do
    end do
    
    do imo = no1+1,no2
       do jmo = no1,imo-1
          do amo = nv1+1,nv2
             do bmo = nv1,amo-1
                kcsf = kcsf + 1
                strb(imo,kcsf) = amo
                strb(jmo,kcsf) = bmo
             end do
          end do
       end do
    end do
    
    !--------- a1/b1 ----------
1000 do imo = no1,no2
       do amo = nv1,nv2
          do jmo = no1,no2         
             do bmo = nv1,nv2
                kcsf = kcsf + 1
                strb(imo,kcsf) = amo
                stra(jmo,kcsf) = bmo
             end do
          end do
       end do
    end do
    
    !--------- a1/b0 .and. a0/b1 ----------    
    do imo = no1,no2
       do amo = nv1,nv2
          kcsf = kcsf + 1
          stra(imo,kcsf) = amo
       end do
    end do
    
    do imo = no1,no2
       do amo = nv1,nv2
          kcsf = kcsf + 1
          strb(imo,kcsf) = amo
       end do
    end do
        
    do kcsf = 1,ncsf
       call Bubble_Sort(stra(:,kcsf),no2,ipa)
       call Bubble_Sort(strb(:,kcsf),no2,ipb)
    end do
    
end subroutine CI_State_Close
!=======================================================================
subroutine Num_Csf_Close(nwo,nwv,ncsf,mcsf)

    implicit none
    integer nwo,nwv,ncsf,mcsf
    
    ncsf = 1 + 2*nwo*nwv + nwo*(nwo-1)*nwv*(nwv-1)/2 + &
            nwo*nwo*nwv*nwv                                 ! CISD
    mcsf = 1 + nwo*(nwo-1)*nwv*(nwv-1)/2 + nwo*nwo*nwv*nwv  ! CID
    
end subroutine Num_Csf_Close
!=======================================================================
!=======================================================================
subroutine CIS_Open(nstate,SCov,EnA,EnB,Ee,g_MoA,g_MoB,g_MoAB,EI,MoCuA, &
           MoCuB,Dipole,N,NoccA,NoccB,E_Nuc)

    implicit none 
    integer nstate,SCov,ii,N,NoccA,NoccB,ncis,EI(N,N)
    integer i,j,a,b,imo,amo,kcis,icis,jcis
    integer nwoA,nwvA,no1A,no2A,nv1A,nv2A
    integer nwoB,nwvB,no1B,no2B,nv1B,nv2B
    integer,allocatable :: str1(:),str2(:)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) PtransA(N,N),PtransB(N,N),sqrDtrans,MoCuA(N,N),MoCuB(N,N),Dipole(N,N,3)
    real(kind=8) EnA(N),EnB(N),Ee,E_Nuc,Ptrans(N,N)
    real(kind=8),allocatable :: CIS_Ham(:,:),Ecis(:),Vcis(:,:),Dtrans(:,:),osc_str(:)
    character(len=13),allocatable :: Ccis(:)
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
       
    ncis = nwoA*nwvA + nwoB*nwvB
    write(*,*)
    write(*,*) " Ncis =", ncis
    
    write(100,*)"************************************************************"
    write(100,*)"Excited states from <AA,BB:AA,BB> singles matrix:           "
    write(100,*)"************************************************************"  
    write(100,*)
    
    allocate(CIS_Ham(ncis,ncis),Ecis(ncis),Vcis(ncis,ncis),str1(ncis), &
             str2(ncis),Ccis(ncis),Dtrans(3,nstate),osc_str(nstate))
             
    !=======================Construct CIS Matrix========================
    kcis = 0
    do imo = no1A,no2A
       do amo = nv1A,nv2A
          kcis = kcis + 1
          str1(kcis) = imo
          str2(kcis) = amo
          write(CTemp1,"(i4)") imo
          write(CTemp2,"(i4)") amo
          Ccis(kcis) = CTemp1//" --> "//CTemp2
       end do
    end do
    
    do imo = no1B,no2B
       do amo = nv1B,nv2B
          kcis = kcis + 1
          str1(kcis) = imo
          str2(kcis) = amo
          write(CTemp1,"(i4)") imo
          write(CTemp2,"(i4)") amo
          Ccis(kcis) = CTemp1//" --> "//CTemp2
       end do
    end do
    !---------------------------alpha-alpha-----------------------------
    do icis = 1,nwoA*nwvA
        i = str1(icis)
        a = str2(icis)
        do jcis = icis,nwoA*nwvA      
            j = str1(jcis)
            b = str2(jcis)
            
            CIS_Ham(icis,jcis) = g_MoA(EI(i,a),EI(j,b)) - g_MoA(EI(i,j),EI(a,b))
            if( i==j .and. a==b ) CIS_Ham(icis,jcis) = CIS_Ham(icis,jcis) + EnA(a) - EnA(i)
            
        end do
    end do
    !----------------------------beta-beta------------------------------
    do icis = nwoA*nwvA+1,ncis
        i = str1(icis)
        a = str2(icis)
        do jcis = icis,ncis     
            j = str1(jcis)
            b = str2(jcis)
            
            CIS_Ham(icis,jcis) = g_MoB(EI(i,a),EI(j,b)) - g_MoB(EI(i,j),EI(a,b))
            if( i==j .and. a==b ) CIS_Ham(icis,jcis) = CIS_Ham(icis,jcis) + EnB(a) - EnB(i)
            
        end do
    end do
    !----------------------------alpha-beta-----------------------------
    do icis = 1,nwoA*nwvA
        i = str1(icis)
        a = str2(icis)
        do jcis = nwoA*nwvA+1,ncis
            j = str1(jcis)
            b = str2(jcis)
            
            CIS_Ham(icis,jcis) = g_MoAB(EI(j,b),EI(i,a))
            
        end do
    end do
    !-------------------------------------------------------------------
    do icis = 1,ncis
        do jcis = icis,ncis
            CIS_Ham(jcis,icis) = CIS_Ham(icis,jcis)
        end do
    end do
    !===================================================================

    call SY_Eigen_System(CIS_Ham,Ecis,Vcis,ncis)
    
    !--------------------- oscillator strength -------------------------
    do ii = 1,nstate
        PtransA = 0.0d0
        PtransB = 0.0d0
        
        do icis = 1,nwoA*nwvA
            i = str1(icis)
            j = str2(icis)
            PtransA =  PtransA + Vcis(icis,ii)*matmul(MoCuA(:,i:i),transpose(MoCuA(:,j:j)))
        end do
        
        do icis = nwoA*nwvA + 1,ncis
            i = str1(icis)
            j = str2(icis)
            PtransB =  PtransB + Vcis(icis,ii)*matmul(MoCuB(:,i:i),transpose(MoCuB(:,j:j)))
        end do
        
        Ptrans = PtransA + PtransB
        Dtrans(1,ii) = sum(Ptrans*Dipole(:,:,1))
        Dtrans(2,ii) = sum(Ptrans*Dipole(:,:,2))
        Dtrans(3,ii) = sum(Ptrans*Dipole(:,:,3))
        sqrDtrans=sum(Dtrans(:,ii)**2)
        osc_str(ii)=2.0d0/3.0d0*Ecis(ii)*sqrDtrans
    end do
    !-------------------------------------------------------------------
    do ii = 1,nstate
    
        write(100,"(A14,I4,A2,F10.4,A3,A2,F10.2,A3,A2,A4,F8.4)") " Excited State",ii," ,",Ecis(ii)*27.2114d0, &
                    " ev"," ,",45.5633d0/Ecis(ii)," nm"," ,"," f= ",osc_str(ii)
        write(100,*) "    Alpha State:"
        do i = 1,nwoA*nwvA
             if( dabs(Vcis(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") Ccis(i),Vcis(i,ii)
        end do
    
        write(100,*) "    Beta State:"
        do i = nwoA*nwvA + 1,ncis
             if( dabs(Vcis(i,ii)) > 10.0d0**(-SCov)) write(100,"(A20,F14.5)") Ccis(i),Vcis(i,ii)
        end do
        write(100,"(A10,F15.8)") "E(cis) = ",Ecis(ii) + Ee + E_Nuc
        write(100,*)  
        
    end do
         
    write(100,*)"************************************************************"      
      

end subroutine CIS_Open
!=======================================================================
subroutine CID_Open(f_MoA,f_MoB,g_MoA,g_MoB,g_MoAB,EI,N,NoccA,NoccB,E_Nuc,Ecid)

    implicit none 
    integer N,NoccA,NoccB,ncsf,mcsf,EI(N,N)
    integer nwoA,nwvA,no1A,no2A,nv1A,nv2A
    integer nwoB,nwvB,no1B,no2B,nv1B,nv2B
    integer n_mem,n_numb
    integer,allocatable :: stra(:,:),strb(:,:)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) Ecid,f_MoA(N,N),f_MoB(N,N),E_Nuc
    real(kind=8),allocatable :: CID_Ham(:,:)
    
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
       
    call Num_Csf_Open(nwoA,nwvA,nwoB,nwvB,ncsf,mcsf)
    call Read_Mem(n_mem,n_numb)
    
    if ( sqrt(real(n_numb)) < real(mcsf) ) then
        write(*,*)
        write(*,*) " Ncid =", mcsf 
        write(*,*)
        write(*,*) " Error : insufficient virtual memory for cid!"
        write(100,*) " Error : insufficient virtual memory for cid!"
        Ecid = - E_Nuc
        return
    end if
    
    write(*,*)
    write(*,*) " Ncid =", mcsf 
    
    allocate(CID_Ham(mcsf,mcsf),stra(no2A,ncsf), &
             strb(no2B,ncsf))
    call CI_State_Open(stra,strb,ncsf,no1A,no2A,nv1A,nv2A,no1B,no2B,nv1B,nv2B)
    call CISD_Hamilton_Open(N,CID_Ham,mcsf,no2A,no2B,stra,strb,g_MoA, &
                              g_MoB,g_MoAB,EI,f_MoA,f_MoB)
      
    call Dig_CI_Ham(CID_Ham,Ecid,mcsf,E_Nuc)

end subroutine CID_Open
!=======================================================================
subroutine CISD_Open(f_MoA,f_MoB,g_MoA,g_MoB,g_MoAB,EI,N,NoccA,NoccB,E_Nuc,Ecisd)

    implicit none 
    integer N,NoccA,NoccB,ncsf,mcsf,EI(N,N)
    integer nwoA,nwvA,no1A,no2A,nv1A,nv2A
    integer nwoB,nwvB,no1B,no2B,nv1B,nv2B
    integer n_mem,n_numb
    integer,allocatable :: stra(:,:),strb(:,:)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2),g_MoAB(N*(N+1)/2,N*(N+1)/2)
    real(kind=8) Ecisd,f_MoA(N,N),f_MoB(N,N),E_Nuc
    real(kind=8),allocatable :: CISD_Ham(:,:)
    
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
        
    call Num_Csf_Open(nwoA,nwvA,nwoB,nwvB,ncsf,mcsf) 
    call Read_Mem(n_mem,n_numb)
    
    if ( sqrt(real(n_numb)) < real(mcsf) ) then
        write(*,*)
        write(*,*) " Ncisd =", ncsf 
        write(*,*)
        write(*,*) " Error : insufficient virtual memory for cisd!"
        write(100,*) " Error : insufficient virtual memory for cisd!"
        Ecisd = - E_Nuc
        return
    end if
    
    write(*,*)
    write(*,*) " Ncisd =", ncsf  
    
    allocate(CISD_Ham(ncsf,ncsf),stra(no2A,ncsf), &
             strb(no2B,ncsf))
    call CI_State_Open(stra,strb,ncsf,no1A,no2A,nv1A,nv2A,no1B,no2B,nv1B,nv2B)
    call CISD_Hamilton_Open(N,CISD_Ham,ncsf,no2A,no2B,stra,strb,g_MoA, &
                              g_MoB,g_MoAB,EI,f_MoA,f_MoB) 
         
    call Dig_CI_Ham(CISD_Ham,Ecisd,ncsf,E_Nuc)

end subroutine CISD_Open
!=======================================================================
subroutine CISD_Hamilton_Open(N,CI_Ham,ncsf,no2A,no2B,stra,strb,g_MoA, &
                              g_MoB,g_MoAB,EI,f_MoA,f_MoB)

    implicit none
    integer N,i,j,kcsf1,kcsf2,ncsf,no2A,no2B,ndif,ndifa,ndifb,ipa,ipb
    integer iobt(no2A+no2B),p,q,r,s,EI(N,N)
    integer stra(no2A,ncsf),strb(no2B,ncsf)
    integer strxa(no2A),strya(no2A),strma(no2A),strna(no2A)
    integer strxb(no2B),stryb(no2B),strmb(no2B),strnb(no2B)
    real(kind=8) g_MoA(N*(N+1)/2,N*(N+1)/2),g_MoB(N*(N+1)/2,N*(N+1)/2) 
    real(kind=8) g_MoAB(N*(N+1)/2,N*(N+1)/2),f_MoA(N,N),f_MoB(N,N)
    real(kind=8) CI_Ham(ncsf,ncsf)
    real(kind=8) Temp1,Temp2,Temp3,Temp4
    
    strxa = 0 ; strya = 0 ; strma = 0 ; strna = 0
    strxb = 0 ; stryb = 0 ; strmb = 0 ; strnb = 0
    CI_Ham = 0.0d0
    
    do kcsf1 = 1,ncsf
       do kcsf2 = kcsf1,ncsf
       
           call Match(no2A,stra(:,kcsf1),stra(:,kcsf2),strxa,strya, &
                      strma,strna,ndifa,ipa)
           call Match(no2B,strb(:,kcsf1),strb(:,kcsf2),strxb,stryb, &
                      strmb,strnb,ndifb,ipb)
           ndif = ndifa + ndifb

           if ( ndif > 2 ) cycle
           !------------------------------------------------------------
           if ( ndif == 2 ) then

               if ( ndifa == 2 ) then
                   p = strxa(1) ; q = strxa(2)
                   r = strya(1) ; s = strya(2)
                   CI_Ham(kcsf2,kcsf1) = g_MoA(EI(s,q),EI(r,p))-g_MoA(EI(r,q),EI(s,p))
                   
               end if
               
               if ( ndifb == 2 ) then
                   p = strxb(1) ; q = strxb(2)
                   r = stryb(1) ; s = stryb(2)
                   CI_Ham(kcsf2,kcsf1) = g_MoB(EI(s,q),EI(r,p))-g_MoB(EI(r,q),EI(s,p))
                   
               end if
               
               if ( ndifa == 1 ) then
                   p = strxa(1) ; q = strxb(1)
                   r = strya(1) ; s = stryb(1)
                   CI_Ham(kcsf2,kcsf1) = g_MoAB(EI(s,q),EI(r,p))
                                     
               end if
               CI_Ham(kcsf2,kcsf1) = CI_Ham(kcsf2,kcsf1)*ipa*ipb
           end if
           !------------------------------------------------------------
           if ( ndif == 1 ) then
               if( kcsf1 == 1 .or. kcsf2 == 1 ) cycle
               
               if( strxa(1) == 0 ) then   
                   iobt(1:no2A) = stra(:,kcsf1) 
                   iobt(no2A+1:no2A+no2B) = strb(:,kcsf1) 
                   p = strxb(1);q = stryb(1)
                   
                   Temp1 = 0.0d0
                   do j = 1,no2A+no2B
                                    
                        i = iobt(j) 
                        if( j <= no2A ) then
                            Temp2 = g_MoAB(EI(p,q),EI(i,i))
                        else
                            Temp2 = g_MoB(EI(i,i),EI(p,q))-g_MoB(EI(q,i),EI(i,p))
                        end if
                         
                        Temp1 = Temp1 + Temp2
                   end do

                   CI_Ham(kcsf2,kcsf1) = f_MoB(q,p) + Temp1

               end if
               
               if( strxb(1) == 0 ) then   
                   iobt(1:no2A) = stra(:,kcsf1) 
                   iobt(no2A+1:no2A+no2B) = strb(:,kcsf1) 
                   p = strxa(1);q = strya(1)
                   
                   Temp1 = 0.0d0
                   do j = 1,no2A+no2B 
                                    
                        i = iobt(j) 
                        if( j > no2A ) then
                            Temp2 = g_MoAB(EI(i,i),EI(p,q))
                        else
                            Temp2 = g_MoA(EI(i,i),EI(p,q))-g_MoA(EI(q,i),EI(i,p))
                        end if
                         
                        Temp1 = Temp1 + Temp2
                   end do

                   CI_Ham(kcsf2,kcsf1) = f_MoA(q,p) + Temp1

               end if
               
               CI_Ham(kcsf2,kcsf1) = CI_Ham(kcsf2,kcsf1)*ipa*ipb
           end if
           !------------------------------------------------------------
           if ( ndif == 0 ) then
           
                iobt(1:no2A) = stra(:,kcsf1) ; iobt(no2A+1:no2A+no2B) = strb(:,kcsf1)
               
                Temp1 = 0.0d0
                do i = 1,no2A
                     p = iobt(i)
                     Temp1 = Temp1 + f_MoA(p,p)
                end do
                
                do i = no2A+1,no2A+no2B
                     p = iobt(i)
                     Temp1 = Temp1 + f_MoB(p,p)
                end do
               
                Temp2 = 0.0d0
                do i = 1,no2A+no2B  
                    p = iobt(i)
                    Temp3 = 0.0d0
                    do j = 1,no2A+no2B                   
                         q = iobt(j) 
                         if( i == j ) cycle 

                         if( j <= no2A .and. i > no2A ) then
                             Temp4 = 0.5d0*g_MoAB(EI(p,p),EI(q,q))
                         else if( i <= no2A .and. j > no2A ) then
                             Temp4 = 0.5d0*g_MoAB(EI(q,q),EI(p,p)) 
                         else if( i <= no2A .and. j <= no2A ) then
                             Temp4 = 0.5d0*(g_MoA(EI(q,q),EI(p,p))-g_MoA(EI(q,p),EI(p,q)))
                         else
                             Temp4 = 0.5d0*(g_MoB(EI(q,q),EI(p,p))-g_MoB(EI(q,p),EI(p,q)))
                         end if
                         
                         Temp3 = Temp3 + Temp4
                    end do
                    Temp2 = Temp2 + Temp3
                end do
               
                CI_Ham(kcsf2,kcsf1) = Temp1 + Temp2
           end if   
           CI_Ham(kcsf1,kcsf2) = CI_Ham(kcsf2,kcsf1)    
       end do
    end do
    
    write(*,*) " CI Hamilton complete."
                   
end subroutine CISD_Hamilton_Open
!=======================================================================
subroutine CI_State_Open(stra,strb,ncsf,no1A,no2A,nv1A,nv2A,no1B,no2B,nv1B,nv2B)

    implicit none
    integer no1A,no2A,nv1A,nv2A,no1B,no2B,nv1B,nv2B
    integer imo,jmo,amo,bmo,ncsf,kcsf
    integer ipa,ipb
    integer stra(no2A,ncsf),strb(no2B,ncsf)

    !--------Ground state--------
    do imo = 1,no2A
       stra(imo,:) = imo
    end do

    do imo = 1,no2B
       strb(imo,:) = imo
    end do
    
    kcsf = 1 
    
    if(no1A == no2A .or. nv1A == nv2A ) goto 999
    !--------- a2/b0 .and. a0/b2 ----------
    do imo = no1A+1,no2A
       do jmo = no1A,imo-1
          do amo = nv1A+1,nv2A
             do bmo = nv1A,amo-1
                kcsf = kcsf + 1
                stra(imo,kcsf) = amo
                stra(jmo,kcsf) = bmo
             end do
          end do
       end do
    end do
    
999 if(no1B == no2B .or. nv1B == nv2B ) goto 1000
    
    do imo = no1B+1,no2B
       do jmo = no1B,imo-1
          do amo = nv1B+1,nv2B
             do bmo = nv1B,amo-1
                kcsf = kcsf + 1
                strb(imo,kcsf) = amo
                strb(jmo,kcsf) = bmo
             end do
          end do
       end do
    end do
    
    !--------- a1/b1 ----------
1000 do imo = no1A,no2A
       do amo = nv1A,nv2A
          do jmo = no1B,no2B         
             do bmo = nv1B,nv2B
                kcsf = kcsf + 1
                stra(imo,kcsf) = amo
                strb(jmo,kcsf) = bmo
             end do
          end do
       end do
    end do
    
    !--------- a1/b0 .and. a0/b1 ----------    
    do imo = no1A,no2A
       do amo = nv1A,nv2A
          kcsf = kcsf + 1
          stra(imo,kcsf) = amo
       end do
    end do
    
    do imo = no1B,no2B
       do amo = nv1B,nv2B
          kcsf = kcsf + 1
          strb(imo,kcsf) = amo
       end do
    end do
        
    do kcsf = 1,ncsf
       call Bubble_Sort(stra(:,kcsf),no2A,ipa)
       call Bubble_Sort(strb(:,kcsf),no2B,ipb)
    end do
    
end subroutine CI_State_Open
!=======================================================================
subroutine Num_Csf_Open(nwoA,nwvA,nwoB,nwvB,ncsf,mcsf)

    implicit none
    integer nwoA,nwvA,nwoB,nwvB,ncsf,mcsf
    
    ncsf = 1 + nwoA*nwvA + nwoB*nwvB + nwoA*(nwoA-1)*nwvA*(nwvA-1)/4 + &
            nwoB*(nwoB-1)*nwvB*(nwvB-1)/4 + nwoA*nwoB*nwvA*nwvB   ! CISD
    mcsf = 1 + nwoA*(nwoA-1)*nwvA*(nwvA-1)/4 + &
            nwoB*(nwoB-1)*nwvB*(nwvB-1)/4 + nwoA*nwoB*nwvA*nwvB   ! CID
    
end subroutine Num_Csf_Open
!=======================================================================

