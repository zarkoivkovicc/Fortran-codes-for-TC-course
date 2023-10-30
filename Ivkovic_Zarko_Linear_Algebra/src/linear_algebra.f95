module linear_algebra
    use inout
    use kinds
    contains
    subroutine Gauss_Jordan(input,output)
        real(rkind),dimension(:,:), intent(in) :: input
        real(rkind),dimension(:,:), intent(out) :: output
        output = input
        call gauss(output)
        call jordan(output)
    end subroutine Gauss_Jordan
    subroutine gauss(matrx)
        !> Gaussian elimination for matrix of a linear system of N equations and N unkowns
        real(rkind),dimension(:,:) :: matrx
        integer :: N
        integer :: i,j ! counters
        N = int(size(matrx,dim=1,kind=rkind))
        call check_gauss(matrx,1)
        matrx(1,:) = matrx(1,:)/matrx(1,1)
        do i=1,N-1
            do j=i+1,N
                matrx(j,:) = matrx(j,:) - matrx(j,i)*matrx(i,:)
                if(j==i+1) then
                matrx(j,:) = matrx(j,:) * 1.0/matrx(j,i+1)
                end if
            end do
            call check_gauss(matrx,i+1)
        end do
    end subroutine gauss
    subroutine jordan(matrx)
        real(rkind),dimension(:,:) :: matrx
        integer :: N
        integer :: i,j ! counters
        N = int(size(matrx,dim=1,kind=rkind))
        do i=N,2,-1
            do j=i-1,1,-1
                matrx(j,:) = matrx(j,:) - matrx(j,i)*matrx(i,:)
            end do
        end do
    end subroutine jordan
    subroutine check_gauss(matri,col)
        !A step in gaussian algorithm, checks if rows should be swapped
        real(rkind),dimension(:,:) :: matri
        real(rkind),dimension(:),allocatable :: temp
        integer :: col, i
        allocate(temp(size(matri,dim=2,kind=rkind)))
        i = col
        do while(matri(i,col).eq. 0)
            i = i + 1
        end do
        temp = matri(i,:); matri(i,:) = matri(col,:); matri(col,:) = matri(i,:)
    end subroutine check_gauss
    subroutine gen_sys(points,system,m)
        !> Generate system of linear equations in order to fit data to m-order polynomial
        implicit none
        real(rkind),dimension(:,:), intent(in) :: points
        real(rkind),dimension(:,:), intent(out) :: system
        real(rkind) :: sum1, sum2
        integer, intent(in) :: m
        integer :: i=1,j=1,k=1 ! counters
        integer :: n ! number of points
        n = int(size(points,dim=1,kind=rkind))
        system(1,1) = n
        do i=1,m+1
            do j=1,i
                sum1 = 0; sum2 = 0
                do k=1,n
                    sum1 = sum1 + points(k,1)**(i-1)
                    sum2 = sum2 + points(k,1)**(2*m-i+1)
                end do
                system(j,i+1-j) = sum1
                system(m+2-j,m+1-i+j) = sum2
            end do
            sum1=0
            do k=1,n
                sum1 = sum1 + points(k,1)**(i-1) * points(k,2)
            end do
            system(i,m+2) = sum1
        end do
        end subroutine gen_sys
    subroutine solve(solution,coef)
            ! Extracts the coefficients from solution matrix and writes them in a file coeff.data
            real(rkind),dimension(:,:), intent(in) :: solution
            real(rkind),dimension(:), intent(out) :: coef
            integer :: m !degree of polynomial
            integer :: i=1,j=1 ! counter
            m = int(size(solution,dim=1,kind=rkind)) - 1
            do i =1,m+1
                do while(solution(i,j) .eq. 0.0)
                    j = j + 1
                end do
                coef(j) = solution(i,m+2)
                j = 1
            end do
            call ofile('coeff.data')
            do i=1,m+1
                write(21,"(i2,f10.5)") i-1, coef(i)
            end do
            close(21)
    end subroutine solve
    function polynom(coef,x) result(y)
        ! Returns P(x) where coef is an array of coeffitients
        real(rkind), dimension(:), intent(in) :: coef
        real(rkind), intent(in) :: x
        real(rkind) :: y
        integer :: m, i !m degree, i is a counter
        m = int(size(coef,dim=1,kind=rkind)) - 1
        y = 0.0
        do i=0,m
            y = y + coef(i+1) * x**i
        end do
    end function polynom
    function r_2(points,coef)
        ! Calculates R^2 coeffitient for data represented by points and polynomial with coeff coef
        real(rkind),dimension(:,:), intent(in) :: points
        real(rkind),dimension(:), intent(in) :: coef
        real(rkind) :: sres, stot, mean, r_2
        integer :: m,n,i !Polynomial of degree m, n points, i is a counter
        m = int(size(coef,dim=1,kind=rkind)) - 1
        n = int(size(points,dim=1,kind=rkind))
        sres = 0; stot = 0; mean = 0
        mean  = sum(points(:,2)) / n
        do i=1,n
            sres = sres + (points(i,2) - polynom(coef,points(i,1)))**2.0
            stot = stot + (points(i,2) - mean)**2.0
        end do
        r_2 = 1.0 - sres/stot
    end function r_2
    subroutine write_pol(coef,points,step)
        !writes a large number of points in the interval to plot the fitted polynomial
        real(rkind),dimension(:,:), intent(in) :: points
        real(rkind),dimension(:), intent(in) :: coef
        real(rkind) :: step, a ,b
        integer :: i, n_step
        a = minval(points(:,1))
        b = maxval(points(:,1))
        n_step = int(2/step + (b - a)/step)
        call ofile('fit.data')
        do i=1,n_step
            write(21,*) a - 1.0 + i*step, polynom(coef,a - 1.0 + i*step)
        end do
        close(21)
    end subroutine write_pol
end module linear_algebra
