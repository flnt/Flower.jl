
"""
Butler-Volmer model, relation between exchande current and potential at electrode, with effects of concentration
c0_H2: reference concentration

```math
i=i_0\\left[\\frac{c_{H_2}}{c_{H_2,0}}^{1/2}\\frac{c_{KOH}}{c_{KOH,0}}\\exp{\\left(\\frac{\\alpha_aF\\eta}{R_\\mu T}\\right)}
-\\frac{c_{H_2O}}{c_{H_2O,0}}\\exp{\\left(\\frac{\\alpha_c F\\eta}{R_\\mu T}\\right)}\\right]
```

```math
\\eta = \\phi_{\text{wall}} - \\phi_{\text{liquid}}
```

"""
function butler_volmer_concentration(alpha_a,alpha_c,c_H2,c0_H2,c_H2O,c0_H2O,c_KOH,c0_KOH,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    # eta = phi_ele1 - phi_ele
    # i_current = i0*(sqrt(c_H2/c0_H2)*(c_KOH/c0_KOH)*exp(alpha_a*Faraday*eta/(Ru*temperature0))-(c_H2O/c0_H2O)*exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    i_current = i0*(sqrt(c_H2/c0_H2)*(c_KOH/c0_KOH)*exp(alpha_a*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0))-(c_H2O/c0_H2O)*exp(-alpha_c*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0)))
    return i_current
end


"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration

alpha_a: transfer coefficients

```math
i=i_0\\left[\\exp{\\left(\\frac{\\alpha_aF\\eta}{R_\\mu T}\\right)}
-\\exp{\\left(\\frac{\\alpha_c F\\eta}{R_\\mu T}\\right)}\\right]
```

"""
function butler_volmer_no_concentration(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    # eta = phi_ele1 - phi_ele
    # i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    i_current = i0*(exp(alpha_a*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0))-exp(-alpha_c*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0)))
    return i_current
end


"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration

alpha_a: transfer coefficients

```math
i=i_0\\left[\\exp{\\left(\\frac{\\alpha_aF\\eta}{R_\\mu T}\\right)}
-\\exp{\\left(\\frac{\\alpha_c F\\eta}{R_\\mu T}\\right)}\\right]
```

"""
function derivative_butler_volmer_no_concentration(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0)
    # eta = phi_ele1 - phi_ele
    # i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    i_butler_derivative = i0*(-alpha_a*Faraday/(Ru*temperature0) * exp(alpha_a*Faraday * (phi_ele1 - phi_ele)/(Ru*temperature0)) - alpha_c*Faraday/(Ru*temperature0) * exp(-alpha_c*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0)) )
    return i_butler_derivative
end

"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration

# Arguments
- `alpha_a::Float64`: transfer coefficient
- `alpha_c::Float64`: transfer coefficient

"""
function butler_volmer_no_concentration!(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0,i_current)
    # eta = phi_ele1 - phi_ele
    # i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    i_current = i0*(exp(alpha_a*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0))-exp(-alpha_c*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0)))
end


"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration

# Arguments
- `alpha_a::Float64`: transfer coefficient
- `alpha_c::Float64`: transfer coefficient

"""
function butler_volmer_no_concentration_concentration_Neumann!(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0,diffusion_coeff,inv_stoechiometric_coeff,a0)
    # eta = phi_ele1 - phi_ele
    # i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    a0 = inv_stoechiometric_coeff*i0*(exp(alpha_a*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0))-exp(-alpha_c*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0)))/(Faraday*diffusion_coeff)
end


"""
Butler-Volmer model, relation between exchande current and potential at electrode without concentration

# Arguments
- `alpha_a::Float64`: transfer coefficient
- `alpha_c::Float64`: transfer coefficient

"""
function butler_volmer_no_concentration_concentration_Neumann(alpha_a,alpha_c,Faraday,i0,phi_ele,phi_ele1,Ru,temperature0,diffusion_coeff,inv_stoechiometric_coeff)
    # eta = phi_ele1 - phi_ele
    # i_current = i0*(exp(alpha_a*Faraday*eta/(Ru*temperature0))-exp(-alpha_c*Faraday*eta/(Ru*temperature0)))
    return inv_stoechiometric_coeff*i0*(exp(alpha_a*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0))-exp(-alpha_c*Faraday*(phi_ele1 - phi_ele)/(Ru*temperature0)))/(Faraday*diffusion_coeff)
end

"""
electrical_conductivity!

Returns electrical conductivity, based on concentration number 2
"""
function electrical_conductivity!(num,concentration,temperature,elec_cond)
    elec_cond .= 2*num.Faraday^2 .*concentration.*num.diffusion_coeff[2]./(num.Ru.*temperature) 
end




"""
butler_volmer_no_concentration_potential_Neumann

Returns Neumann boundary condition for electrical potential
"""
function butler_volmer_no_concentration_potential_Neumann!(num,phi_eleD,concentration,temperature,a0)
                   
    # a0 = num.i0*(exp(num.alpha_a*num.Faraday*(num.phi_ele1 - phi_eleD)/(num.Ru*temperature)) - exp(-num.alpha_c*num.Faraday*(num.phi_ele1 - phi_eleD)/(num.Ru*temperature)))/(2*num.Faraday^2 *concentration*num.diffusion_coeff[2]/(num.Ru*temperature) )
    
    a0 = butler_volmer_no_concentration_potential_Neumann_no_struct!(num.i0, num.alpha_a, num.alpha_c, num.Faraday, num.phi_ele1, phi_eleD, 
    num.Ru, num.diffusion_coeff[2], concentration, temperature)
end

"""
butler_volmer_no_concentration_potential_Neumann

Returns Neumann boundary condition for electrical potential
"""
function butler_volmer_no_concentration_potential_Neumann(num,phi_eleD,concentration,temperature)
                   
    # a0 = num.i0*(exp(num.alpha_a*num.Faraday*(num.phi_ele1 - phi_eleD)/(num.Ru*temperature)) - exp(-num.alpha_c*num.Faraday*(num.phi_ele1 - phi_eleD)/(num.Ru*temperature)))/(2*num.Faraday^2 *concentration*num.diffusion_coeff[2]/(num.Ru*temperature) )
    
    return butler_volmer_no_concentration_potential_Neumann_no_struct!(num.i0, num.alpha_a, num.alpha_c, num.Faraday, num.phi_ele1, phi_eleD, 
    num.Ru, num.diffusion_coeff[2], concentration, temperature)
end


# function butler_volmer_no_concentration_potential_Neumann(num,phi_eleD,concentration,temperature,a0)
                   
#     a0 .= num.i0*(exp.(num.alpha_a*num.Faraday.*(num.phi_ele1 .- phi_eleD)./(num.Ru.*temperature)) .- exp.(-num.alpha_c.*num.Faraday.*(num.phi_ele1 .- phi_eleD)./(num.Ru.*temperature)))./(2*num.Faraday^2 .*concentration.*num.diffusion_coeff[2]./(num.Ru.*temperature) )
# end

function butler_volmer_no_concentration_potential_Neumann_no_struct!(i0,alpha_a,alpha_c,Faraday,phi_ele1,phi_eleD,Ru,diffusion_coeff,concentration,temperature)
                   
    a0 = i0*(exp(alpha_a*Faraday*(phi_ele1 - phi_eleD)/(Ru*temperature)) - exp(-alpha_c*Faraday*(phi_ele1 - phi_eleD)/(Ru*temperature)))/(2*Faraday^2 *concentration*diffusion_coeff/(Ru*temperature))

end


"""
Computes conductivity based on a concentration of electrolyte 
"""
function compute_ele_cond(Faraday,diffusion_coeff,Ru,temperature,concentration)
   
    return 2*Faraday^2 *concentration*diffusion_coeff/(Ru*temperature)
end

