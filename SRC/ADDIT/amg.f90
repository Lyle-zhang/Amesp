program hmg

    implicit none
    integer narg 
    integer N,i,j,k,l,Charge,Spin_Mult,rec_char
    real(kind=8),allocatable:: Rn(:,:)
    character(len=80) Ch0,Ch1,Atoms,method,option,outs,bas_set
    character(len=200) fname,sfname
    character(len=2),allocatable:: Atom(:)
    
    narg = IARGC()
    if( narg > 0 ) then
       call getarg(1,fname)
    else
       write(*,*) " input file name,or drag the file into terminal directly. "
       read(*,*) fname
    end if
    
    call Filename(fname,sfname,rec_char)
    open(unit = 10,file = fname)   
    read(10,"(A1)") Ch0
    close(10)

    if ( Ch0 == "#" .or. Ch0 == "%" ) then
    
        open(unit = 20,file = fname)
        open(unit=100,file = sfname(1:rec_char)//".aip")
        
        k = -1
        do while(.true.)
           read(20,"(A1)") Ch1
           k = k + 1
           if( Ch1 == "#" ) exit
        end do
        
        read(20,*)
        read(20,*)
        read(20,*)
        read(20,*)
        
        i = 0
        do while(.true.)
           read(20,"(A2)") Atoms
           i = i + 1
           if( Atoms == "  " ) exit
        end do
        N = i - 1
        allocate(Rn(3,N),Atom(N))
        close(20) 
        
        open(unit = 30,file = fname)
        do i = 1,k
            read(30,*)
        end do
        
        read(30,"(A2,A30)") Ch1,method
        read(30,*)
        read(30,*)
        read(30,*)
        read(30,*) Charge,Spin_Mult
        do i = 1,N
           read(30,*) Atom(i),Rn(1:3,i)
        end do
        close(30)
        
        do i = 1,20
           if( method(i:i) == "/" ) exit
        end do
        j = i

        do i = j+1,20
           if( method(i:i) == " " ) exit
        end do
        l = i
        
        outs = " out=1"
        write(100,"(G0,I3,I3)") N,Charge,Spin_Mult
        write(100,"(A7,A12,A6,A6)") method(j+1:l),method(1:j-1)," none",outs  
        write(100,"(A41)") "maxcyc= 128 conver= 8 diis=on guess=atden"   
        write(100,*)   
        do i = 1,N
             write(100,"(A3,F20.8,2F14.8)") Atom(i),Rn(1:3,i)
        end do
        
        write(*,*)
        write(*,*) " Transform gaussian input file to amesp input file successfully."
    end if
    !-------------------------------------------------------------------
    if( Ch0 /= "#" .and. Ch0 /= "%" ) then
    
        open(unit = 20,file = fname)
        open(unit=100,file = sfname(1:rec_char)//".com")
        read(20,*) N,Charge,Spin_Mult
        allocate(Rn(3,N),Atom(N))
        read(20,*) bas_set,method,option,outs
        read(20,*)
        read(20,*)
        do i = 1,N
           read(20,*) Atom(i),Rn(1:3,i)
        end do
        close(20)

        do i = 1,20
           if( method(i:i) == " " ) exit
        end do
        j = i    

        do i = 1,20
           if( bas_set(i:i) == " " ) exit
        end do
        k = i 
                    
        write(100,"(A3,A20,A18)") "#p ",method(1:j-1)//"/"//bas_set(1:k-1), &
                                  " pop=full nosymm" !gfinput iop(3/33=1) A37 
        write(100,*)
        write(100,"(A10)") "Warm_Cloud"
        write(100,*)
        write(100,"(G0,I3)") Charge,Spin_Mult
        
        do i = 1,N
            write(100,"(A3,F20.8,2F14.8)") Atom(i),Rn(1:3,i)
        end do
        write(100,*)
        write(100,*)
        
        write(*,*)
        write(*,*) " Transform amesp input file to gaussian input file successfully."        
    end if
    close(100)
    
end program hmg
!=======================================================================
subroutine Filename(fname,sfname,rec_char)

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
