!Zarko Ivkovic, University of Barcelona
program Huckel
    use kinds;use inout, only: read_xyz, ofile; use math, only: build_H,diagonalize,sort_nor_eig, fill_orbitals, mulliken
    implicit none
    real(rkind), dimension(:,:), allocatable :: H, eigenvec, bond ! Hamiltonian, Eigenvectors matrix and bond order matrix
    real(rkind), dimension(:), allocatable :: eigenval, occ, pi_occ ! eigenvalues, occupation numbers and pi occupation numbers
    integer, dimension(:,:), allocatable :: amatrix ! adjacency matrix
    real(rkind) :: alpha, beta ! huckel parameters
    integer :: charge, n, i,j ! charge, number of carbon atoms, counters
    character(20) :: name ! name of the system
    read(5,*) name, charge
    call read_xyz(amatrix,trim(name)//'.xyz') !Read coordinates
    alpha = -11.4; beta = -0.8;
    n = int(size(amatrix,dim=1,kind=rkind))
    allocate(H(n,n)) ; allocate(eigenvec(n,n)) ; allocate(eigenval(n)); allocate(occ(n))
    call build_H(amatrix,H,alpha,beta) ! Build Huckel hamiltonian based on the coordinates
    ! Calculate everything
    call diagonalize(H,eigenvec,eigenval)
    call sort_nor_eig(eigenvec,eigenval)
    call fill_orbitals(eigenval,occ,charge)
    allocate(pi_occ(n)); allocate(bond(n,n))
    call mulliken(eigenvec,occ,bond,pi_occ,amatrix)
    !PRODUCE THE OUTPUT
    call ofile(trim(name)//'.evec')
    write(21,*) "EIGENVECTORS"
    do i=1,n
        write(21,*) (eigenvec(i,j), j=1,n)
    end do
    close(21)
    call ofile(trim(name)//'.eval')
    write(21,*) "EIGENVALUES"
    write(21,*) (eigenval(i), i=1,n)
    close(21)
    call ofile(trim(name)//'.occ')
    write(21,*) "Occupation numbers"
    write(21,*) (occ(i), i=1,n)
    close(21)
    call ofile(trim(name)//'.piocc')
    write(21,*) " PI Occupation numbers"
    write(21,*) (pi_occ(i), i=1,n)
    close(21)
    call ofile(trim(name)//'.bo')
    write(21,*) "Bond orders"
    do i=1,n
        write(21,*) (bond(i,j), j=1,n)
    end do
    close(21)
end program Huckel
