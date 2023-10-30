!Zarko Ivkovic, University of Barcelona
program Newton_Raphson
    use kinds; use optimization, only: grad, hess, inverse
    implicit none
    real(rkind), parameter :: Conv = 1E-8 ! Convergence criteria
    real(rkind), dimension(2) :: gradient ! gradient
    real(rkind), dimension(2) :: coord,old_coord ! coordinate and coordinate from the previous step
    
    real(rkind), dimension(2,2) :: hessian, inv ! hessian matrix and it's inverse
    integer :: p = 1 ! counter
    real(rkind) :: length, diff
    coord(1) = 1.0; coord(2) = 3.0 ! initial guess
    length =  norm2(grad(coord(:),func)) !error as the norm of gradient
    diff = 1 ! making sure we enter the loop
    write(6,*) "       x         y           f        grad(x)      grad(y)&
    &      H(xx)       H(xy)       H(yx)       H(yy)      Error"
    do while (length > Conv .and. diff > Conv ) !Convergence is reached either when a) gradient is below a treshold or b) the value of function between two succesive points differs less than a treshold
        gradient = grad(coord(:),func)
        hessian = hess(coord,func)
        inv = inverse(hessian)
        length = norm2(grad(coord(:),func)) !norm of the gradient
        old_coord = coord !save the previous point
        coord = coord - matmul(inv,gradient) !calculate the next point
        diff = abs(func(coord(1),coord(2)) - func(old_coord(1),old_coord(2)))
        write(6,"(10(es12.4))") old_coord(1), old_coord(2), func(old_coord(1),old_coord(2)),gradient(1)&
        ,gradient(2), hessian(1,1), hessian(1,2), hessian(2,1), hessian(2,2), diff
        p = p + 1
    end do
    contains
    function func(x,y) result(f)
        !> Returns the f(x) = sin(x+y) + (x-y)^2 - 1.5x + 3.5y + 3
        real(rkind), intent(in) :: x,y
        real(rkind) :: f
        f = dsin(x+y) + (x-y)**2.0 - 1.5*x + 3.5*y + 3.0
    end function func
end program Newton_Raphson