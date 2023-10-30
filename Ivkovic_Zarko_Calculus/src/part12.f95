! Zarko Ivkovic, University of Barcelona
program Romberg
use kinds; use integrals, only : Romb
implicit none
real(rkind),dimension(10,10) :: R !Romberg matrix 10x10
real(rkind) :: a = 1.0 ! lower limit of integration
real(rkind) :: b = 3.0 ! upper limit of integration
real(rkind) :: Result ! to store the final converged value
real(rkind) :: conv = real(1E-8,rkind) ! convergence criteria
logical :: check = .true. ! is convergence criteria NOT REACHED
integer :: j,k ! counters
integer :: t=1,y=1 ! Position of converged result in the field
R(1,1) = (b-a)/2*(func(a)+func(b)) ! Initialize the first element
write(6,"(f13.9)") R(1,1) ! print the first element
do k=2,10 ! run over rows
    do j=1,k ! run over columns
        if(j == 1) then ! if it's the first row, calculate elemement using helper function romb
            R(k,j) = 0.5*R(k-1,j) + Romb(k,func,a,b) 
        else ! if not, calculate it using recursive relation from prevous elements
            R(k,j) = R(k,j-1) + (R(k,j-1) - R(k-1,j-1))/(4**(j-1)-1)
        end if
        if(abs(R(k,j)-R(k,j-1)) < conv .and. check .eqv. .true.) then ! check for convergence
            t = k; y = j; Result = R(k,j) ! store the position and value of the element that fulfills the convergnce criteria
            check=.false.
        end if
    end do
    write(6,"(55(f13.9))") (R(k,j),j=1,k) ! write the row of triangular matrix
end do
write(6,"(a,i1,a,i1,a,f13.8)") "The converged element is R(",t,",",y,") = ", Result ! Print the result
contains
function func(x) result(retval)
    !> Function that's integrated
    !> In this case, it's sin(x^2) - cos(2x)
    use kinds
    real(rkind), intent(in) :: x
    real(rkind) :: retval
    retval = sin(x**2) - cos(2*x)
end function func
end program Romberg