program main
        use globals
        use ziggurat
        use MD

        integer :: seed,i
        logical :: es
        real(kind=8), allocatable :: Es_min(:)
        real(kind=8) :: Ek, inst_T, p_in, p

        inquire(file='seed.dat',exist=es)
        if(es) then
            open(unit=10,file='seed.dat',status='old')
            read(10,*) seed
            close(10)
            print *,"  * Leyendo semilla de archivo seed.dat"
        else
            seed = 24583490
        end if

        call zigset(seed)
       
        
        rho = 0.2*sigma**(-3)
        L = 10.0
        N = rho*(L**3)
        print *, N
        r_cut = sigma*2.5
        v_cut = 4*eps*(-(sigma/r_cut)**6+(sigma/r_cut)**12)
        dt = 0.01
        T = 3.0


        allocate(r(3,N))
        allocate(v(3,N))
        allocate(F(3,N))

        call init_coords()
        p_in = update_E_and_F()
        call update_lgv_F()
        
        Es_min = E_minimization(5000,1)
        call initiate_velocities()
        
        open(unit=10, file = "produccion.dat", action="write", status="replace", form="formatted")
        write(10,"(A34)") "step,E total,Ep,Ek,F max,v max,T,p"

        do i=1,10000
                p = integrate()
                if(mod(i,10)==0) then
                        Ek = get_Ek()
                        inst_T = get_T()
                        write(10,"(I5,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8,A1,f16.8)") &
                        i,",",E_tot+Ek,",",E_tot,",",Ek,",",maxval(F),",",maxval(v),",",inst_T,",",p
                        write(*,"(A6,I5)") "step: ", i
                        print *, E_tot + Ek 
                        call save_coords("produccion.xyz")
                end if
        end do
        close(10)
                
end program main







