!Zarko Ivkovic, University of Barcelona
! NOTE: I am not using kinds module here because the given subroutine uses real*8 and integers
program Gau_Leg
    use integrals, only : Gauss_int
    implicit NONE
    real(8), dimension(1:2) :: Quad ! Values of two succesive quadratures
    real(8), parameter :: lower = 1.0, upper = 3.0, Conv = 1D-8 !lower limit, upper limit, convergence criteria
    real(8) :: diff ! difference between two succesive steps
    integer :: Nsub,p !number of intervals, number of steps
    p=1; Nsub=2
    Quad(2) = Gauss_int(lower,upper,Nsub,func)
    write(6,"(a,i2,a,i2,a,es15.8,a)") "Step: ",p, " Number of quadrature points: " ,Nsub,&
    "  Integral: ",Quad(2),"  Difference: "
    diff = 1 ! making sure that we enter the loop
    do while (abs(diff) > Conv .and. (Nsub < 12)) ! Unit the convergence is reached or the number of points is larger than 12
        Nsub = Nsub + 1 !Increase the number of points
        Quad(1) = Quad(2) ! Save the previous step
        Quad(2) = Gauss_int(lower,upper,Nsub,func)
        diff = Quad(2) - Quad(1) !Calculate difference between two succesive quadratures
        p = p+1
        write(6,"(a,i2,a,i2,a,es15.8,a,es11.4)") "Step: ",p, " Number of quadrature points: " ,Nsub,&
    "  Integral: ",Quad(2),"  Difference: ",diff
    end do
    write(6,*) "FINAL OUTPUT"
    write(6,"(a,i2)") "Number of quadrature points: ", Nsub
    write(6,"(a,es15.8)") "Integral: ",Quad(2)
    contains
    function func(x) result(retval)
        !> Function that's integrated
        !> In this case, it's sin(x^2) - cos(2x)
        real(8), intent(in) :: x
        real(8) :: retval
        retval = dsin(x**2) - dcos(2*x)
    end function func
end program Gau_Leg