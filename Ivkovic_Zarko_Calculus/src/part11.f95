! Zarko Ivkovic, University of Barcelona
program simpson_integration
    use kinds; use integrals, only: simpson
    implicit none
    real(rkind), dimension(1:2) :: I ! Values of two succesive quadratures
    real(rkind), parameter :: a = 1.0, b = 3.0, Eps = 1E-8 !lower limit, upper limit, convergence criteria
    real(rkind) :: diff ! difference between two succesive steps
    integer(ikind) :: N,j !number of intervals, number of steps
    j=1; N=1
    I(2) = simpson(a,b,N,func)
    write(6,"(a,i4,a,i6,a,es15.8,a)") "Step: ",j, " Number of subintervals: " ,N,&
    "  Integral: ",I(2),"  Difference: "
    diff = 1 ! making sure that we enter the loop
    do while (abs(diff) > Eps) ! Unit the convergence is reached
        N = N*2 !double the number of subintervals
        I(1) = I(2)
        I(2) = simpson(a,b,N,func)
        diff = I(2) - I(1)
        j = j+1
        write(6,"(a,i4,a,i6,a,es15.8,a,es11.4)") "Step: ",j, " Number of subintervals: " ,N,&
    "  Integral: ",I(2),"  Difference: ",diff
    end do
    write(6,*) "FINAL OUTPUT"
    write(6,"(a,i4)") "Number of interations needed: ",j
    write(6,"(a,i10)") "Number of subintervals: ",N
    write(6,"(a,i8)") "Number of abcissa points: ", N+1
    write(6,"(a,es15.8)") "Integral: ",I(2)
    contains
        function func(x) result(retval)
            !> Function that's integrated
            !> In this case, it's sin(x^2) - cos(2x)
            use kinds
            real(rkind), intent(in) :: x
            real(rkind) :: retval
            retval = sin(x**2) - cos(2*x) 
        end function func
end program simpson_integration