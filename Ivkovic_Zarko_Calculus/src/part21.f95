program Steepest_Descent
    use kinds; use optimization, only : grad
    implicit none
    real(rkind), parameter :: step = 0.3, Conv = 1E-8 ! Step in the optimization algorithm and conv criteria
    real(rkind), dimension(2) :: gradient, n_gradient ! gradient and normalized gradient
    integer(ikind),parameter :: maxiter = 40 ! maximal number of iterations
    integer(ikind) :: i = 1 ! counter
    real(rkind),dimension(maxiter+1,2) :: coord ! coordinates of points visited during the optimization
    real(rkind) :: diff !norm of the distance between two succesive points
    coord(1,1) = 1; coord(1,2) = 3 ! initial guess
    diff =  1 ! making sure we enter the loop
    write(6,*) "       x         y            f       grad(x)      grad(y)    ngrad(x)     ngrad(y)    Error"
    do while ((i <= maxiter) .and. (diff > Conv .or. norm2(grad(coord(i,:),func)) > Conv) ) !Converge criteria
        gradient = grad(coord(i,:),func)
        n_gradient = gradient / norm2(gradient) !Normalize gradient
        write(6,"(8(es12.3))") coord(i,1), coord(i,2), func(coord(i,1),coord(i,2)),gradient(1)&
        ,gradient(2),n_gradient(1),n_gradient(2), diff
        coord(i+1,1) = coord(i,1) - n_gradient(1)*step !make a step
        coord(i+1,2) = coord(i,2) - n_gradient(2)*step ! make a step
        diff = func(coord(i+1,1),coord(i+1,2)) - func(coord(i,1),coord(i,2)) ! error as the difference between the function at two succesive points
        i = i + 1
    end do
    call write_plot(coord,'data')
    contains
    function func(x,y) result(f)
        !> Returns the f(x) = sin(x+y) + (x-y)^2 - 1.5x + 3.5y + 3
        real(rkind), intent(in) :: x,y
        real(rkind) :: f
        f = sin(x+y) + (x-y)**2.0 - 1.5*x + 3.5*y + 3.0
    end function func
    subroutine write_plot(A,file)
        ! DESCRIPTION: 
        !>  Writes a matrix A to a file named <file>.
        !>  The resulting file will be structured like this:
        !>  x1 y1 f(x1,y1)
        !>  x2 y2 f(x2,y2)
        !>        . . .
        !>  xn yn f(xn,yn)
        !------------------------------------------------------------------------------
        implicit none
        real(rkind), dimension(:,:) :: A
        character(*) :: file
        integer :: n
        call ofile(file)
        n = size(A,dim=1)
        do i=1,n
            write(21,*) A(i,1), A(i,2), func(A(i,1),A(i,2))
        end do
        close(21)
    end subroutine write_plot
    subroutine ofile(name)
        ! DESCRIPTION: 
        !< Opens a new file on unit 21. If a file with the same name already exists, it gets replaced.
        implicit none
        character(*), intent(in) :: name
        integer :: ierr
      
        open (21, FILE=name, STATUS='unknown', iostat=ierr)
        if ( ierr == 0) then
            close(21, status='delete')
            open(21, FILE=name, STATUS='new')
        endif
      end subroutine ofile
end program Steepest_Descent