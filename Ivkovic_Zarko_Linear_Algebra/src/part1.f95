! Zarko Ivkovic, University of Barcelona
program Linear_System
    use kinds ; use linear_algebra; use inout
    implicit none
    real(rkind), dimension(:,:), allocatable :: system, solution, points
    real(rkind), dimension(:), allocatable :: polyn
    integer :: m,n !We want to fit n point to m degree polinomial
    character(3) :: set ! Number of set of points
    read(5,*) set, m
    call read_points(points,'points_'//trim(set)//'.data')
    n = int(size(points,dim=1,kind=rkind))
    allocate(system(m+1,m+2)) ; allocate(solution(m+1,m+2)); allocate(polyn(m+1))
    call gen_sys(points,system,m)
    call Gauss_Jordan(system,solution)
    call solve(solution,polyn)
    write(6,*) set, r_2(points,polyn)
    call write_pol(polyn,points,real(0.1,rkind))
end program Linear_System
