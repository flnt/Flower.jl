"""
Computes average value at interface for scalar
"""
function mean_intfc_non_null(scalD,iscal,grid,iLS)

    index = iLS+1
    num=0
    nonzero = 0.0

    # cf veci @view a[g.ny*g.nx*(p-1)+1:g.ny*g.nx*p]

    for i in grid.ny*grid.nx*(index-1)+1:grid.ny*grid.nx*index
        if abs(scalD[i,iscal]) .> 0.0
            nonzero += scalD[i,iscal]
            num += 1
        end
    end

    if num == 0
        print("\n no intfc in mean_intfc_non_null")
        return 0
    else
        nonzero /= num
        return nonzero
    end
    
end