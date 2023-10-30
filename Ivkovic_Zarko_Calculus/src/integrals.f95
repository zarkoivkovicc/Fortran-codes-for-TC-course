module integrals
    use kinds
    contains
    function simpson(c,d,Nit,func) result(retval)
        !> Calculates integral of func from a to b (a<b) using Simpson's rule with Nit subintervals
        use kinds
        implicit none
        interface
        function func(x)
            use kinds
            real(rkind), intent(in) :: x
            real(rkind) :: func
        end function func
        end interface
        
        
        real(rkind), intent(in) :: c,d ! lower and upper limits of integration
        real(rkind) :: retval ! quadrature
        real(rkind) :: sum_e, sum_o ! sums with even and odd members that appear in the formula
        integer(ikind), intent(in) :: Nit ! Number of subintervals
        integer(ikind) :: k ! Simple counter
        real(rkind) :: h ! Step size
        sum_e = 0; sum_o = 0 ! Initialize sums to 0
        h = (d-c)/(2*Nit) ! Calculate step size
        do k=2,(2*Nit-2),2 !sum of even elements
            sum_e = sum_e + func(c + k*h)
        end do
        do k=1,(2*Nit-1),2 ! sum of odd elements
            sum_o = sum_o + func(c + k*h)
        end do
        retval = h*(func(c) + func(d) + 2*sum_e + 4*sum_o)/3.0 ! simpson composite rule
    end function simpson
    function Romb(n,func,c,d) result(res)
        !> Helper function to calculate the elements of the first column
        use kinds
        interface
        function func(x)
            use kinds
            real(rkind), intent(in) :: x
            real(rkind) :: func
        end function func
        end interface
        real(rkind) :: res,h ! result, step of integration
        real(rkind), intent(in) :: c,d ! integration limints
        integer(ikind), intent(in) :: n ! Row that we are calculating the first element of
        integer(ikind) :: i !counter
        h = (d-c)/(2**(n-1)) ! integration step
        do i=1,(2**(n-2)) ! Do the integration
            res = res + func(c + (2*i - 1)*h)
        end do
        res = res*h
    end function Romb
    function Gauss_int(a,b,Nin,func) result(retval)
        !> Calculates integral of func from a to b (a<b) using Gauss-Legendere interpolation
        implicit none
        interface
        function func(x)
            use kinds
            real(8), intent(in) :: x
            real(8) :: func
        end function func
        end interface
        
        real(8), intent(in) :: a,b ! lower and upper limits of integration
        real(8) :: retval ! quadrature
        real(8) :: c,m
        integer, intent(in) :: Nin ! Number of points
        integer :: k ! Simple counter
        real(8), dimension(12) :: weight, abcisa ! arrays of weights and abcisa points
        weight = 0; abcisa = 0 ! Initialize arrays to 0
        c = (a+b)*0.5D0 !calculate coffeicents from the formula
        m = (b-a)*0.5D0
        call sub_GauLeg(real(-1.0,8),real(1.0,8),abcisa,weight,Nin) !Calculate weights and quadrature points in the interval [-1,1]
        retval = 0
        do k=1,Nin
            retval = retval + weight(k)*func(c + m*abcisa(k)) ! The function with change of variables
        end do
        retval = retval * m
    end function Gauss_int
    SUBROUTINE sub_GauLeg(X1,X2,t,w,n)
        !beste modu bat X1 eta X2 kendu, eta XL eta XM
        !bukaera t(i) = -z eta t(n+1-i)= +z izango dira
    
          IMPLICIT NONE
          
          INTEGER :: m,i,j
          INTEGER,intent(in) :: n   !Number of Gaussian points
          REAL(8), dimension(n),intent(out) :: W,t
          REAL(8) :: XM,XL,X1,X2,EPS,P1,P2,P3,pi,Z1,Z,PP
    
        !Relative precision
           EPS = 1.D-14
    
          !double precision arccosine. Pi value=3.14159
         pi = DACOS(-1.D0)
    
        !N = number of Gauss Points
        !Roots are symmetric in the interval - so only need to find half of them  
           m = (n + 1) / 2
        
        !The coats are going to be  -1 and 1, Gauss-Legendre 
        !Variable change
              XM=0.5D0*(X1+X2)
              XL=0.5D0*(X2-X1)
    
    
        !Loop over the desired roots
              DO i = 1,m
             Z = DCOS (pi * (i - 0.25D0)/(n + 0.5D0))
        !Starting with the above approximation to the i-th root,
        !we enter the main loop of refinement by NEWTON'S method   
     10      P1 = 1.D0
             P2 = 0.D0
    
        !Loop up the recurrence relation to get the Legendre
        !polynomial evaluated at z                
             DO j = 1,n
                P3 = P2
                P2 = P1
                P1 = ((2.D0 * j - 1.D0) * Z * P2 - (j - 1.D0) * P3)/j
             END DO
    !p1 is now the desired Legendre polynomial.
    !We next compute pp, its derivative, by a standard relation involving also p2, 
    !the polynomial of one lower order. 
             PP = n * (Z * P1 - P2)/(Z * Z - 1.D0)
             Z1 = Z
             Z = Z1 - P1/PP    ! Newton's Method  */
    
             IF (DABS(Z-Z1) .GT. EPS) GO TO 10
    
        ! Roots will be symmetric about the origin  
             t(i) = XM - XL * Z
             t(n + 1 - i) = XM + XL * Z
        !Compute the weight and its symmetric counterpart 
             W(i) = 2.D0 * XL/((1.D0 - Z * Z) * PP * PP)
             W(n + 1 - i) = W(i)
          END DO  
    
        RETURN   !not neccesary
          END SUBROUTINE Sub_GauLeg
end module integrals
