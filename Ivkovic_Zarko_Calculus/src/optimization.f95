module optimization
    use kinds
    contains
    function grad(point,fun)
        ! Takes a point in 2-D space and returns gradient and normalized gradient of a function fun in that point  
        ! point: array of two numbers (x,y)
        ! returns (gradx,grady)
        interface
        function fun(x,y) result(f)
            use kinds
            real(rkind), intent(in) :: x,y
            real(rkind) :: f
        end function fun
        end interface
        !> Calculate gradient of f(x,y) at point, where point is two dimensional array (x,y)
        real(rkind), dimension(2), intent(in) :: point !array of coordinates of a point: (x,y)
        real(rkind), dimension(2) :: grad !gradient vector (array of 2 numbers: grad_x and grad_y)
        grad(1) = der_x(point,fun)
        grad(2) = der_y(point,fun)
    end function grad
    function hess(point,fun)
        ! Takes a point in 2-D space and returns gradient and normalized gradient of a function fun in that point  
        ! point: array of two numbers (x,y)
        ! returns (gradx,grady)
        interface
        function fun(x,y) result(f)
            use kinds
            real(rkind), intent(in) :: x,y
            real(rkind) :: f
        end function fun
        end interface
        !> Calculate gradient of f(x,y) at point, where point is two dimensional array (x,y)
        real(rkind), dimension(2), intent(in) :: point !array of coordinates of a point: (x,y)
        real(rkind) :: h !step size for numerical calculation of the gradient 
        real(rkind), dimension(2,2) :: hess !gradient vector (array of 2 numbers: grad_x and grad_y)
        real(rkind), dimension(3,2) :: steps
        h = real(1E-6,rkind) ! OPTIMAL STEP SIZE
        steps(1,1) = h  ; steps(1,2) = 0 ! step for second derivative with respect to x
        steps(2,1) = 0  ; steps(2,2) = h !step for second derivative with respect to y
        hess(1,1) = (der_x(point + steps(1,:),fun) - der_x(point - steps(1,:),fun))/(2.0*h)
        hess(2,2) = (der_y(point + steps(2,:),fun) - der_y(point - steps(2,:),fun))/(2.0*h)
        hess(2,1) = (der_x(point + steps(2,:),fun) - der_x(point - steps(2,:),fun))/(2.0*h)
        hess(1,2) = hess(2,1) ! Hessian matrix is symmetric
    end function hess
    function der_x(point,fun)
        ! Calculates partial derivative with respect to x of a function f(x,y) in point  
        interface
        function fun(x,y) result(f)
            use kinds
            real(rkind), intent(in) :: x,y
            real(rkind) :: f
        end function fun
        end interface
        real(rkind), dimension(2), intent(in) :: point !array of coordinates of a point: (x,y)
        real(rkind) :: h, der_x !step size for numerical calculation of the derivative, derivative
        h = real(1E-6,rkind) ! OPTIMAL STEP SIZE
        der_x = (fun(point(1)+h,point(2))-fun(point(1)-h,point(2)))/(2.0*h)
    end function der_x
    function der_y(point,fun)
        ! Calculates partial derivative with respect to x of a function f(x,y) in point  
        interface
        function fun(x,y) result(f)
            use kinds
            real(rkind), intent(in) :: x,y
            real(rkind) :: f
        end function fun
        end interface
        real(rkind), dimension(2), intent(in) :: point !array of coordinates of a point: (x,y)
        real(rkind) :: h, der_y !step size for numerical calculation of the derivative, derivative
        h = real(1E-6,rkind)
        der_y = (fun(point(1),point(2)+h)-fun(point(1),point(2)-h))/(2.0*h)
    end function der_y
    function inverse(matrx)
        ! Calculates inverse of a function using Gauss-Jordan method
        real(rkind),dimension(:,:) :: matrx
        real(rkind),dimension(:,:),allocatable :: inverse
        real(rkind),dimension(:,:),allocatable :: expded ! expanded matrix
        integer :: N
        integer :: i,j,k ! counters
        N = int(size(matrx,dim=1,kind=rkind))
        allocate(expded(N,2*N))
        allocate(inverse(N,N))
        inverse = 0.0
        expded = 0.0
        expded(1:N,1:N) = matrx; 
        do k=1,N !form expanded matrix
            expded(k,k+N) = 1.0
        end do
        !Gaussian elimination
        expded(1,:) = expded(1,:) / expded(1,1) ! make the first element 1
        do i=1,N-1
            do j=i+1,N
                expded(j,:) = expded(j,:) - expded(j,i)*expded(i,:) !make 0
                if(j==i+1) then
                expded(j,:) = expded(j,:) * 1.0/expded(j,i+1) !if it's diagonal make it 1
                end if
            end do
        end do
        !Back-substitution
        do i=N,2,-1
            do j=i-1,1,-1
                expded(j,:) = expded(j,:) - expded(j,i)*expded(i,:)
            end do
        end do
        inverse = expded(:,N+1:2*N)
        deallocate(expded)
    end function inverse
end module optimization