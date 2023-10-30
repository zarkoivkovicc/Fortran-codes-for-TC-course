module math
    use kinds
    contains
    subroutine build_H(amatrix,H,alpha,beta)
        integer, dimension(:,:), intent(in) :: amatrix !adjacency matrix
        real(rkind), dimension(:,:), intent(out) :: H ! Hamiltonian
        real(rkind), intent(in) :: alpha, beta ! huckel parameters
        integer :: n, i ,j !number of C atoms, counters
        n = int(size(amatrix,dim=1))
        do i=1,n
            H(i,i) = alpha
            do j=i+1,n
                if (amatrix(i,j) == 1) then
                    H(i,j) = beta
                    H(j,i) = beta
                else
                    H(i,j) = 0.0
                    H(j,i) = 0.0
                end if
            end do
        end do
    end subroutine build_H
    subroutine diagonalize(matrix,eigenvec,eigenval)
        implicit none
        real(rkind),dimension(:,:), intent(in) :: matrix ! matrix we are diagonazing
        real(rkind),dimension(:,:), intent(out) :: eigenvec ! matrix in which columns are eigenvectors
        real(rkind),dimension(:), intent(out) :: eigenval ! array of corresponding eigenvalues
        real(rkind), dimension(:,:), allocatable :: O, temp !O is unitary matrix we are using as a step in jacobi algorithm, temp is temporary matrix where we do all tranformations
        real(rkind) :: atg ! theta in formulas
        integer,dimension(2) :: pivot !the element we are setting to 0
        integer :: n, i
        n = int(size(matrix,dim=1,kind=rkind))
        allocate(O(n,n)); allocate(temp(n,n))
        call indentity(O) ; call indentity(eigenvec)
        temp = matrix
        pivot = maxloc(abs(ignore_diagonal(matrix))) !find the largest off-diagonal elements
        do while(abs(temp(pivot(1),pivot(2))) >= 1E-12)
            if (abs(temp(pivot(1),pivot(1)) - temp(pivot(2),pivot(2)))>= 1E-14) then !check if the theta is pi/4 or calculate it
                atg = 0.5*ATAN(2.0*temp(pivot(1),pivot(2))/(temp(pivot(1),pivot(1)) - temp(pivot(2),pivot(2))))
            else 
                atg = 0.5*ATAN(1D0)
            end if
            O(pivot(1),pivot(1)) = COS(atg)
            O(pivot(2),pivot(1)) = SIN(atg)
            O(pivot(1),pivot(2)) = -O(pivot(2),pivot(1))
            O(pivot(2),pivot(2)) = O(pivot(1),pivot(1))
            temp = matmul(matmul(transpose(O),temp),O)
            eigenvec = matmul(eigenvec,O)
            pivot = maxloc(abs(ignore_diagonal(temp)))
            call indentity(O)
        end do
        do i=1,n
            eigenval(i) = temp(i,i)
        end do
        deallocate(temp); deallocate(O)
    end subroutine diagonalize
    subroutine indentity(matrix)
        !> Makes an identity matrix out of matrix
        implicit none
        real(rkind),dimension(:,:), intent(inout) :: matrix
        integer :: n,i
        n = int(size(matrix,dim=1,kind=rkind))
        matrix = 0
        do i=1,n
            matrix(i,i) = 1
        end do
    end subroutine indentity
    function ignore_diagonal(matrix) result(res)
        !> Returns the same matrix with all diagonal elements set to 0
        implicit none
        real(rkind),dimension(:,:), intent(in) :: matrix
        real(rkind),dimension(:,:),allocatable :: res
        integer :: n,i
        n = int(size(matrix,dim=1,kind=rkind))
        allocate(res(n,n))
        res = matrix
        do i=1,n
            res(i,i) = 0
        end do
    end function ignore_diagonal
    subroutine sort_nor_eig(eigenvec,eigenval)
        !> Sorts all eigenvectors by increasing eigenvalue, and sorts eigenvalues so the correspondence is the same
        real(rkind),dimension(:,:), intent(inout) :: eigenvec
        real(rkind),dimension(:), intent(inout) :: eigenval
        real(rkind),dimension(:), allocatable :: temp_v
        real(rkind) :: temp
        integer :: i,j,n
        n = int(size(eigenval,dim=1,kind=rkind))
        allocate(temp_v(n))
        do i=1,n
            do j=i+1,n
                if(eigenval(j) < eigenval(i)) then
                    temp = eigenval(i)
                    temp_v = eigenvec(:,i)
                    eigenval(i) = eigenval(j)
                    eigenvec(:,i) = eigenvec(:,j)
                    eigenval(j) = temp
                    eigenvec(:,j) = temp_v
                end if
            end do
        end do
        do i=1,n
            eigenvec(:,i) = eigenvec(:,i) / norm2(eigenvec(:,i))
        end do
    end subroutine sort_nor_eig
    subroutine fill_orbitals(eigenval,occ,charge)
        !> Fills orbitals taking into account AUFBAU, Pauli and Huckel principle
        real(rkind),dimension(:), intent(in) :: eigenval
        real(rkind),dimension(:), intent(out) :: occ
        integer :: charge, n_e, n, i, deg, j
        n = int(size(eigenval,dim=1,kind=rkind))
        n_e = n - charge
        i = 1; deg = 1; occ = 0
        do while(n_e > 0)
            occ(i) = occ(i) + 1
            n_e = n_e - 1
            i = i + 1
            if(abs(eigenval(i)-eigenval(i-1)) > 1E-8 .and. n_e > 0) then
                occ(i-1) = occ(i-1) + 1
                n_e = n_e - 1
            else if(n_e > 0) then
                deg = i-1
                do while(abs(eigenval(i)-eigenval(i-1)) <= 1E-8 .and. n_e > 0)
                    occ(i) = occ(i) + 1
                    n_e = n_e - 1
                    i = i + 1
                end do
                do j=deg,i
                    if(n_e > 0) then
                        occ(j) = occ(j) + 1
                        n_e = n_e - 1
                    end if
                end do
            end if
        end do
    end subroutine fill_orbitals
    subroutine mulliken(eigenvec,occ,bond,pi_occ,amatrix)
        !Performs mulliken analysis
        real(rkind),dimension(:,:), intent(in) :: eigenvec
        real(rkind),dimension(:), intent(in) :: occ
        real(rkind),dimension(:), intent(out) :: pi_occ
        real(rkind),dimension(:,:), intent(out) :: bond
        integer, dimension(:,:), intent(in) :: amatrix
        integer :: i,j,k,n
        n = int(size(eigenvec,dim=1,kind=rkind))
        pi_occ = 0; bond = 0
        do j=1,n
            do i=1,n
                pi_occ(j) = pi_occ(j) + occ(i)*eigenvec(j,i)**2.0
            end do
        end do
        do j=1,n
            do k=j+1,n
                if (amatrix(j,k) == 1) then
                    do i=1,n
                        bond(j,k) = bond(j,k) + occ(i)*eigenvec(j,i)*eigenvec(k,i)
                        bond(k,j) = bond(j,k)
                    end do
                end if
            end do
        end do
    end subroutine mulliken
end module math