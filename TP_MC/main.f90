program MC
        
        T = 300
        beta = 1/(k_b*T)

        A = get_matrix_rand(N) !NxN 1 y -1 asignados aleatoriamente
        A = load_matrix(dir)   !opcionalmente carga la matriz desde un archivo
        
        E_tot = get_total_energy(A)   !energia computada con una matriz entera
        E_med = 0
        E_sqr = 0
        M_tot = get_total_mag(A)      !magnetizacion desde una matriz entera
        M_med = 0
        M_sqr = 0

        do i=1,n_steps
                x = integer(uni()*N)  !elegimos una posicion aleatoria
                y = integer(uni()*N)
                
                dE,s_k = get_delta_E(x,y,A)  !calcula el delta E por cambiar el spin de posicion x,y de M,
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

                if(i%1000) then       !ver como se hace en fortran
                        
                        write(..., E_tot)   ! !se guardan la energia actual y el promedio hasta este paso, etc 
                                            ! !puede ser en el mismo archivo
                        write(..., E_med/i)
                        write(..., M_tot...)

                        save_matrix(dir,A)  !guardar la matriz
                end if
        end do

end program MC

