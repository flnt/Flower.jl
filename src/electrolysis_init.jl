function Poiseuille_fmax(x,v_inlet_max,L0)
    return 4*v_inlet_max*x/L0*(1-x/L0)
end

function Poiseuille_favg(x,v_inlet_moy,L0)
    return 6*v_inlet_moy*x/L0*(1-x/L0)
end


function test_Poiseuille(num,velD,grid_v)

    vel = reshape(vec1(velD,grid_v), grid_v)

    #error = 0.0
    error = vel .- Poiseuille_fmax.(grid_v.x,num.v_inlet,num.L0)
    error_rel = maximum(abs.(error))/num.v_inlet
    error_rel_min = minimum(abs.(error))/num.v_inlet


    scal_error_border = maximum(abs.(vecb_B(velD,grid_v) .- Poiseuille_fmax.(grid_v.x[1,:],num.v_inlet,num.L0)))
    scal_error_border = scal_error_border/num.v_inlet

    printstyled(color=:red, @sprintf "\n Velocity error bottom: %.2e\n" scal_error_border )

    scal_error_border = maximum(abs.(vecb_T(velD,grid_v) .- Poiseuille_fmax.(grid_v.x[1,:],num.v_inlet,num.L0)))
    scal_error_border = scal_error_border/num.v_inlet

    printstyled(color=:red, @sprintf "\n Velocity error top: %.2e\n" scal_error_border )

    printstyled(color=:red, @sprintf "\n Velocity error bulk: %.2e %.2e %.10e\n" error_rel error_rel_min vel[1,1])

end
