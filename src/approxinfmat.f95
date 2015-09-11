subroutine approxinfmat(p, q, psi, imatrix)
    ! Large-sample information matrix for ARMA process
    ! based on 
    !  Box and Jenkins (1970)
    !  and InformationMatrixARMA in FitARMA R package

    implicit none

    integer, intent(in) :: p, q
    integer :: i, j
    double precision, intent(in), dimension(p + q) :: psi
    !double precision, intent(inout) :: imdet
    double precision, dimension(p) :: ptmp
    double precision, dimension(q) :: qtmp
    double precision, dimension(2 * max(p, q),1) :: tmpmat
    double precision, dimension(p + q,p + q) :: imatrix


    external arcov

    if (p > 0) then                
        ptmp = psi(1:p)          
        call arcov(ptmp, ptmp, p, p, tmpmat(1:(2 * p), 1))
        do i = 1, p
            imatrix(i, i:p) = tmpmat(p:(2 * p - i), 1)
        end do          
    end if
    if (q > 0) then         
        qtmp = -psi((p + 1):(p + q))
        call arcov(qtmp, qtmp, q, q, tmpmat(1:(2 * q), 1))
        do i = 1, q
            imatrix(p + i,(p + i):(p + q)) = tmpmat(q:(2 * q - i), 1)
        end do          
    end if        
    if (p > 0 .AND. q > 0) then          
        call arcov(ptmp, qtmp, p, q, tmpmat(1:(p + q), 1))         
        do i = 1, p
            do j = 1, q                           
                imatrix(i, p + j) = -tmpmat(q + i-  j, 1)
            end do
        end do          
    end if        
!    call dpotrf('u', p + q, imatrix, p + q, info)
!    if(info == 0) then !parameter redundancy can cause the information matrix to non-positive definite                  
!        imdet = 1.0d0
!        do i = 1, (p + q)
!            imdet = imdet * imatrix(i, i)
!        end do    
!    end if
end subroutine approxinfmat