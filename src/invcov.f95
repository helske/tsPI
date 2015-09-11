subroutine armainvcov(m,n,zt,tt,rt,p1,invcov)

    implicit none

    integer, intent(in) ::  m, n
    integer ::  t,i

    double precision, intent(inout),dimension(n,n) :: invcov
    double precision, intent(in), dimension(1,m) :: zt  
    double precision, intent(in), dimension(m,m) :: tt 
    double precision, intent(in), dimension(m,1) :: rt
    double precision, intent(in), dimension(m,m) ::  p1

    double precision, dimension(m,m,n+1) :: pt
    double precision, dimension(n) :: ft
    double precision, dimension(m,1,n) :: kt
    double precision, dimension(m,m,n) :: lt
    double precision, dimension(m,m) :: tp
    double precision, dimension(m,1) :: pz

    double precision, dimension(m,1) :: lk,lk2
    double precision, dimension(n,n) :: lowc

    double precision, external :: ddot

    external dgemm
    external dgemv
    external daxpy
    external dsymm
    external dsymv

    kt=0.0d0
    lt=0.0d0
    pt=0.0d0
    ft=0.0d0
    pt(1:m,1:m,1) = p1

    do t = 1, n
        call dsymv('u',m,1.0d0,pt(1:m,1:m,t),m,zt(1,1:m),1,0.0d0,pz,1) ! p symmetric!
        ft(t) = ddot(m,zt(1,1:m),1,pz,1)  !ft 
        lt(1:m,1:m,t) = tt(1:m,1:m)
        call dgemv('n',m,m,1.0d0/ft(t),tt(1:m,1:m),m,pz,1,0.0d0,kt(1:m,1,t),1)
        call dgemm('n','n',m,m,1,-1.0d0,kt(1:m,1,t),m,zt(1,1:m),1,1.0d0,lt(1:m,1:m,t),m) !lt = t - kz
        call dsymm('r','u',m,m,1.0d0,pt(1:m,1:m,t),m,tt(1:m,1:m),m,0.0d0,tp,m)
        call dgemm('n','t',m,m,m,1.0d0,tp,m,lt(1:m,1:m,t),m,0.0d0,pt(1:m,1:m,t+1),m)   
        call dsyr('u',m,1.0d0,rt(1:m,1),1,pt(1:m,1:m,t+1),m)  
    end do


    lowc=0.0d0
    invcov=0.0d0     
    do t = 1, n
        lk(:,1)=kt(:,1,t)     
        do i = (t+1), n      
            lowc(i,t) = -lk(1,1)
            call dgemv('n',m,m,1.0d0,lt(:,:,i),m,lk(:,1),1,0.0d0,lk2(:,1),1)
            lk = lk2             
        end do
        lowc(t,t) = 1.0d0
        invcov(t,t) = 1.0d0/ft(t)
    end do
    lowc(n,n) = 1.0d0
    invcov(n,n) = 1.0d0/ft(n)     
       
    call dtrmm('l','l','t','u',n,n,1.0d0,lowc,n,invcov,n)
    call dtrmm('r','l','n','u',n,n,1.0d0,lowc,n,invcov,n)

end subroutine armainvcov