"""
Computes average value at interface for scalar
"""
function mean_intfc_non_null(scalD,iscal,grid,iLS)

    index = iLS+1
    num=0
    nonzero = 0.0
    ni = grid.nx * grid.ny
    nb = 2 * grid.nx + 2 * grid.ny
    nt = (num.nLS + 1) * ni + nb

    for i in 1:nt
        if abs(veci(scalD[i,iscal],grid,index)) .> 0.0
            nonzero += veci(scalD[i,iscal],grid,index)
            num += 1
        end
    end

    if num == 0
        print("\n no intfc in mean_intfc_non_null")
        return 0
    else
        nonzero ./= num
        return nonzero
    end
    
end