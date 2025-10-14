program MC
        !UNIDADES: relativas, J = kb = 1
        use ising_mc
        use ziggurat
        implicit none
        integer, parameter :: dp = selected_real_kind(15, 307)
        real(kind=dp) :: T, r, E_tot, E_med, E_sqr, M_tot, M_med, M_sqr, dE, beta, acept
        integer :: N, x, y, n_steps,s_k,i
        integer, allocatable :: A(:,:)
        logical :: start

        open(unit=30, file="in.dat", status="old", action="read")
        read(30,*) T,N,n_steps
        close(30)

        beta = 1/T
        

        allocate(A(0:N+1,0:N+1))  !A tiene que ser alocado antes de get_matrix_rand o de load_matrix
        inquire(file="matriz.dat",exist=start)

        if(start) then
                call load_matrix(A,N,'matriz.dat')   !opcionalmente carga la matriz desde un archivo
        else 
                call get_matrix_rand(A,N) !NxN 1 y -1 asignados aleatoriamente
        end if

        E_tot = get_total_energy(A)   !energia computada con una matriz entera
        E_med = 0
        E_sqr = 0
        M_tot = get_total_mag(A)      !magnetizacion desde una matriz entera
        M_med = 0
        M_sqr = 0
        acept = 0
        
        open(unit=30, file="matrix_out.dat", status = "replace", action="write")
        open(unit = 10, file="out.dat", status = "replace", action = "write")
        write(10,*)"Paso,E_tot,E_med,E_sqr,M_tot,M_med,M_sqr,acept_ratio"
        do i=1,n_steps  
                x = int(uni()*N)+1  !elegimos una posicion aleatoria
                y = int(uni()*N)+1
                s_k = A(x,y)
                dE = get_dE(x,y,A,s_k)  !calcula el delta E por cambiar el spin de posicion x,y de M,
                                             !devuelve dE y el valor s_k antes del cambio
                
                if(dE < 0) then
                        A(x,y) = -s_k
                        E_tot = E_tot + dE
                        M_tot = M_tot + (-2*s_k)
                        acept = acept + 1
                else
                        r = uni()
                        if(r < exp(-beta*(dE))) then
                                A(x,y) = -s_k
                                E_tot = E_tot + dE
                                M_tot = M_tot + (-2*s_k)
                                acept = acept + 1
                        end if
                end if
                !actualizamos las filas y columnas de condicion periodica
                A(0,1:N) = A(N,1:N)
                A(N+1,1:N) = A(1,1:N)
                A(1:N,0) = A(1:N,N)
                A(1:N,N+1) = A(1:N,1)

                E_med = E_med + E_tot
                M_med = M_med + M_tot
                E_sqr = E_sqr + (E_tot*E_tot)
                M_sqr = M_sqr + (M_tot*M_tot)
                if (MOD(i,10000) == 0) then       
                        call save_matrix(A,'matriz.dat')  !guardar la matriz
                        write(10,*) i,",",E_tot,",",E_med/i,",",E_sqr/i,",",M_tot,",",M_med/i,",",M_sqr/i,",",acept/i    ! !se guardan la energia actual y el promedio hasta este paso, etc 
                end if
                if(mod(i,10000)==0) then
                        write(30,*) A
                end if
        end do
        close(10)
        close(30)
end program MC

