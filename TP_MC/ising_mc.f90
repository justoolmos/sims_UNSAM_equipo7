module ising_mc

contains

function get_matrix_rand(N) result(matrix)
        use ziggurat
        implicit none
        integer, intent(in) :: N
        integer :: i,j
        integer, allocatable :: matrix(:,:)

        allocate(matrix(0:N+1,0:N+1))
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

        matrix(0,:) = matrix(N,:)
        matrix(N+1,:) = matrix(1,:)
        matrix(:,0) = matrix(:,N)
        matrix(:,N+1) = matrix(:,1)

end function get_matrix_rand


function get_total_energy(A) result(E)
        
        integer, intent(in) :: A(:,:)
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

        integer, intent(in) :: A(:,:)
        real :: M
        integer :: l
        
        l = size(A,1)-1
        M = sum(A(1:l,1:l))

end function get_total_mag

function get_dE(x,y,A,s_k) result(dE)

        integer, intent(in) :: x,y,A(:,:),s_k
        real :: dE

        dE = 2*(s_k*A(x+1,y) + s_k*A(x,y+1) + s_k*A(x-1,y) + s_k*A(x,y-1))
end function get_dE
        

subroutine print_matrix(A)
        integer, intent(in) :: A(:,:)
        integer :: l,i

        l = size(A,1)-2    !N
        do i=1,l
                write(*,"(*(I3))") A(i,:)   !*(I3) significa entero de 3 digitios, el * es para repetirlo cuanto haga falta
        end do
        print *, " "
end subroutine print_matrix


end module ising_mc



