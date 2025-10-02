module ising_mc

contains

subroutine get_matrix_rand(matrix,N) 
        use ziggurat
        implicit none
        integer, intent(in) :: N
        integer :: i,j
        integer, intent(inout) :: matrix(0:,0:)

        matrix = 0
        do j=1,N
                do i=1,N
                        if (uni() < 0.5) then
                                matrix(i,j) = 1
                        else
                                matrix(i,j) = -1

                        end if
                end do
        end do

        matrix(0,1:N) = matrix(N,1:N)
        matrix(N+1,1:N) = matrix(1,1:N)
        matrix(1:N,0) = matrix(1:N,N)
        matrix(1:N,N+1) = matrix(1:N,1)
        
end subroutine get_matrix_rand


function get_total_energy(A) result(E)
        
        integer, intent(in) :: A(0:,0:)
        real :: E
        integer :: i,j,l,s

        l = size(A,1)-2
        E = 0.0

        do j=1,l
                do i=1,l

                        s = A(i,j)
                        E = E + s*A(i-1,j) + s*A(i,j+1)
                end do
        end do
end function get_total_energy


function get_total_mag(A) result(M)

        integer, intent(in) :: A(0:,0:)
        real :: M
        integer :: l
        
        l = size(A,1)-2
        M = sum(A(1:l,1:l))

end function get_total_mag

function get_dE(x,y,A,s_k) result(dE)

        integer, intent(in) :: x,y,A(0:,0:),s_k
        real :: dE

        dE = 2*(s_k*A(x+1,y) + s_k*A(x,y+1) + s_k*A(x-1,y) + s_k*A(x,y-1))
end function get_dE
        

subroutine print_matrix(A)
        integer, intent(in) :: A(0:,0:)
        integer :: l,i

        l = size(A,1)-2    !N
        do i=1,l
                write(*,"(*(I3))") A(i,1:l)   !*(I3) significa entero de 3 digitios, el * es para repetirlo cuanto haga falta
        end do
        print *, " "
end subroutine print_matrix


subroutine print_matrix_full(A)
        integer, intent(in) :: A(0:,0:)
        integer :: l,i

        l = size(A,1)    !N
        do i=0,l-1
                write(*,"(*(I3))") A(i,:)   !*(I3) significa entero de 3 digitios, el * es para repetirlo cuanto haga falta
        end do
        print *, " "
end subroutine print_matrix_full

subroutine save_matrix(A, fname)
  implicit none
  integer, intent(in) :: A(0:,0:)
  character(*), intent(in) :: fname
  integer :: u, i, l
  l = size(A,1)

  open(newunit=u, file=fname, status="replace", action="write", form="formatted")
  do i = 0, l-1
     write(u,"(*(I3))") A(i,:)
  end do
  close(u)
end subroutine save_matrix

subroutine load_matrix(A, n, fname)
  implicit none
  integer, intent(inout) :: A(0:,0:)
  integer, intent(in) :: n
  character(*), intent(in) :: fname
  integer :: u, i

  open(newunit=u, file=fname, status="old", action="read")

  do i = 0 , n+1
     read(u,*) A(i, :)
  end do

  close(u)
end subroutine load_matrix

end module ising_mc



