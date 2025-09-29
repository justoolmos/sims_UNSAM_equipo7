program MC
        !UNIDADES: mks

        implicit none

        real :: T, k_b, E_tot, E_med, E_sqr, M_tot, M_med, M_sqr
        integer :: N, x, y, n_steps
        real, allocatable :: A(:,:)
           
        T = 300
        k_b =  1.3806e-23
        beta = 1/(k_b*T)
        N = 20
        n_steps = 1000

        A = get_matrix_rand(N) !NxN 1 y -1 asignados aleatoriamente
       !A = load_matrix(dir)   !opcionalmente carga la matriz desde un archivo
        
        E_tot = get_total_energy(A)   !energia computada con una matriz entera
        E_med = 0
        E_sqr = 0
        M_tot = get_total_mag(A)      !magnetizacion desde una matriz entera
        M_med = 0
        M_sqr = 0
        
        open(unit = 10, file="out.dat", status = "replace", action = "write")
        write(10,*) "Paso,E_tot,E_med,E_sqr,M_tot,M_med,M_sqr"
        do i=1,n_steps  
                x = integer(uni()*N)  !elegimos una posicion aleatoria
                y = integer(uni()*N)
                
                dE,s_k = get_dE(x,y,A)  !calcula el delta E por cambiar el spin de posicion x,y de M,
                                             !devuelve dE y el valor s_k antes del cambio
                
                if(dE < 0) then
                        A(x:y) = -s_k
                        E_tot = E_tot + dE
                        M_tot = M_tot + (-2*s_k)
                else
                        r = uni()
                        if(r < exp(-beta*(dE))) then
                                A(x:y) = -s_k
                                E_tot = E_tot + dE
                                M_tot = M_tot + (-2*s_k)
                        end if
                end if
                
                E_med = E_med + E_tot
                M_med = M_med + M_tot
                E_sqr = E_sqr + (E_tot*E_tot)
                M_sqr = M_sqr + (M_tot*M_tot)

                if (MOD(i,1000) == 0) then       

                        write(10,*) i,,",",E_tot,",",E_med/i,",",E_sqr/i,",",M_tot,",",M_med/i,",",M_sqr/i    ! !se guardan la energia actual y el promedio hasta este paso, etc 
                                            ! !puede ser en el mismo archivo
                                             
                        save_matrix(dir,A)  !guardar la matriz
                end if
        end do
contains

function get_matrix_rand(N) result(matrix)
        use ziggurat
        implicit none
        integer, intent(in) :: N
        integer :: i,j
        integer, allocatable :: matrix(:,:)

        allocate(matrix(N,N))
        do j=1,N
                do i=1,N
                        if (uni() < 0.5) then
                                matrix(i,j) = 1
                        else
                                matrix(i,j) = -1

                        end if
                end do
        end do

end function get_matrix_rand


function get_total_energy(A) result(E)
        
        integer, intent(in) :: A(:,:)
        real :: E
        integer :: i,j,l,v_1,v_2,s

        l = shape(A)(1)
        E = 0

        do j=1,l
                do i=1,l

                        s = A(i,j)

                        if(i==1) then
                                v_1 = A(l,j)
                        else
                                v_1 = A(i-1,j)
                        end if

                        if(j==l) then
                                v_2 = A(i,1)
                        else
                                v_2 = A(i,j+1)
                        end if

                        E = E + s*v_1 + s*v_2
                end do
        end do
end function get_total_energy


end program MC

