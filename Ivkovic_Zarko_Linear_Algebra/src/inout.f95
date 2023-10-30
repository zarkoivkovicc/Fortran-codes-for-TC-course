module inout
  use kinds
  contains
  !functions is_blank, number_of_columns and number_of_rows are taken from fortran-lang/stdlib repository
function is_blank(c)
    implicit none
    logical :: is_blank
    character(len=1), intent(in) :: c !! The character to test.
    integer :: ic

    ic = iachar(c)             ! TAB
    is_blank = (c == ' ') .or. (ic == int(z'09'));
    return
end function
integer function number_of_columns(s)
    implicit none
    integer,intent(in) :: s
    integer :: ios
    character :: c
    logical :: lastblank

    rewind(s)
    number_of_columns = 0
    lastblank = .true.
    do
      read(s, '(a)', advance='no', iostat=ios) c
      if (ios /= 0) exit
      if (lastblank .and. .not. is_blank(c)) number_of_columns = number_of_columns + 1
      lastblank = is_blank(c)
    end do
    rewind(s)
  end function number_of_columns
  
integer function number_of_rows(s) result(nrows)
  implicit none
  integer,intent(in)::s
  integer :: ios
  real :: r
  complex :: z

  rewind(s)
  nrows = 0
  do
    read(s, *, iostat=ios) r
    if (ios /= 0) exit
    nrows = nrows + 1
  end do
  rewind(s)
  ! If there are no rows of real numbers, it may be that they are complex
  if( nrows == 0) then
    do
      read(s, *, iostat=ios) z
      if (ios /= 0) exit
      nrows = nrows + 1
    end do
    rewind(s)
  end if
  end function number_of_rows

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
subroutine read_points(A,file)
        ! DESCRIPTION: 
        !>  Reads a points from a file and stores it in matrix A.
        !>  File must be structured like this:
        !>  x_1 y_1
        !>  x_2 y_2
        !>        . . .
        !>  x_n y_2
        !------------------------------------------------------------------------------
        implicit none
        real(rkind), dimension(:,:), allocatable :: A
        character(*) :: file
        integer :: i,j,n,m

        open(11, STATUS="OLD", FILE=file)
        n = number_of_rows(11)
        m = number_of_columns(11)
        allocate (A(n,m))
        do i=1,n
            read(11,*) (A(i,j), j=1,m)
        end do
        close(11)
    end subroutine read_points
subroutine read_xyz(amatrix,file)
  ! DESCRIPTION: 
  !>  Reads xyz file to obtain adjacency matrix
  !------------------------------------------------------------------------------
  implicit none
  real(rkind), dimension(:,:), allocatable :: coord, temp !adjacency matrix, coordinates of c atoms and temp matrix
  integer, dimension(:,:), allocatable :: amatrix
  character(*) :: file
  integer :: n,nc = 0 ! total number of atoms and number of carbon atoms
  character(len=1) :: element
  real(rkind) :: dist
  integer :: i,j ! counter
  open(11, STATUS="OLD", FILE=file)
  read(11,*) n
  allocate(temp(n,3)); allocate(coord(n,3))
  read(11,*) ! Title
  do i=1,n
    read(11,*) element, temp(i,1), temp(i,2), temp(i,3)
    if(element == 'C') then ! if elements is Carbon, save it
      nc = nc + 1
      coord(nc,:) = temp(i,:)
    end if
  end do
  close(11)
  deallocate(temp)
  allocate(amatrix(nc,nc))
  do i=1,nc
    do j=i+1,nc
      dist = norm2(coord(i,:)-coord(j,:))
      if (dist >= 1.30 .and. dist <= 1.50) then ! if 2 carbon atoms are around bond distance, count them as bonded
        amatrix(i,j) = 1
      else
        amatrix(i,j) = 0
      end if
    end do
  end do
end subroutine read_xyz
end module inout
