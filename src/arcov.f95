subroutine arcov(ar, ma, p, q, b)
    ! Cross-covariances of auxiliary AR process of ARMA(p, q) model
    ! Based on
    !  McLeod, A.I. (1975),
    !  Derivation of the theoretical autocorrelation function of
    !  autoregressive moving-average time series,
    !  Applied Statistics 24, 255-256.
    !  and tccfAR function in FitARMA R package
    implicit none

    integer, intent(in) ::  p, q
    integer ::  k, i,j, info,imj,ijq
    double precision, intent(in), dimension(p) :: ar
    double precision, intent(in), dimension(q) :: ma
    double precision, intent(out), dimension(p + q, 1) :: b
    double precision, dimension(p + q, p + q) :: x
    integer, dimension(p + q) :: ipiv

    external dgesv


    k = p + q
    x = 0.0d0

    do i = 1, k
        do j = 1, k
            imj = i - j
            ijq = i + j - q -1
            if(i>q) then
                if(i >j .AND. imj .LE.q) then
                    x(i,j) = ma(imj)
                else
                    if(i > q .AND. imj == 0) then
                        x(i,j) = -1.0d0
                    end if
                end if
            else
                if(ijq > 0 .AND. ijq .LE. p) then
                    x(i,j) = ar(ijq)
                else
                    if(ijq==0) then
                        x(i,j) = -1.0d0
                    end if
                end if
            end if
        end do
    end do

    b = 0.0d0
    b(1,1) = -1.0d0
    call dgesv(k, 1, x, k, ipiv, b, k, info)

end subroutine arcov
