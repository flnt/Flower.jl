# Base.show(io::IO, BC::Boundaries) = print(io,"name $(BC.name) \n left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
# Base.show(io::IO, BC::Boundaries) = print(io," left = $(typeof(BC.left)) $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
# Base.show(io::IO, BC::Boundaries) = print(io," left = $(typeof(BC.left)) $(BC.left.val) \n right = $(typeof(BC.right)) $(BC.right.val) \n bottom = $(typeof(BC.bottom)) $(BC.bottom.val) \n top = $(typeof(BC.top)) $(BC.top.val) \n")


# function print_BC_line(BC)
#     # print(" left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
#     print(" left = $(typeof(BC.left)) $(BC.left.val) right = $(typeof(BC.right)) $(BC.right.val) bottom = $(typeof(BC.bottom)) $(BC.bottom.val) top = $(typeof(BC.top)) $(BC.top.val) \n")
# end

# function print_BC(BC)
#     # print(" left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
#     print(" left = $(typeof(BC.left)) $(BC.left.val) \n right = $(typeof(BC.right)) $(BC.right.val) \n bottom = $(typeof(BC.bottom)) $(BC.bottom.val) \n top = $(typeof(BC.top)) $(BC.top.val) \n")
# end


# does not work with precompilation
# Base.show(io::IO, BC::Boundaries) = print(io," left = $(typeofBC(BC.left)) $(BC.left.val) \n right = $(typeofBC(BC.right)) $(BC.right.val) \n bottom = $(typeofBC(BC.bottom)) $(BC.bottom.val) \n top = $(typeofBC(BC.top)) $(BC.top.val) \n")
# Base.show(io::IO, BC::BoundariesInt) = print(io," left = $(typeofBC(BC.left)) $(BC.left.val) \n right = $(typeofBC(BC.right)) $(BC.right.val) \n bottom = $(typeofBC(BC.bottom)) $(BC.bottom.val) \n top = $(typeofBC(BC.top)) $(BC.top.val) \n LS 1 = $(typeofBC(BC.LS[1])) $(BC.LS[1].val) \n")


function print_BC_LS_html(BC,name;io=stdout)
    # ::IO
    # print(" left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
    print(io," <tr> <td> $(name) </td> <td>$(typeofBC(BC.left)) $(BC.left.val)</td> <td>$(typeofBC(BC.right)) $(BC.right.val)</td> <td>$(typeofBC(BC.bottom)) $(BC.bottom.val) </td> <td> $(typeofBC(BC.top)) $(BC.top.val) </td> <td> $(typeofBC(BC.LS[1])) $(BC.LS[1].val) </td> </tr> \n")
end

function print_BC_html(BC,name;io=stdout)
    # ::IO
    # print(" left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
    print(io," <tr> <td> $(name) </td> <td>$(typeofBC(BC.left)) $(BC.left.val)</td> <td>$(typeofBC(BC.right)) $(BC.right.val)</td> <td>$(typeofBC(BC.bottom)) $(BC.bottom.val) </td> <td> $(typeofBC(BC.top)) $(BC.top.val) </td> </tr> \n")
end

function print_BC_line(BC)
    # print(" left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
    print(" left = $(typeofBC(BC.left)) $(BC.left.val) right = $(typeofBC(BC.right)) $(BC.right.val) bottom = $(typeofBC(BC.bottom)) $(BC.bottom.val) top = $(typeofBC(BC.top)) $(BC.top.val) \n")
end

function print_BC(BC)
    # print(" left = $(BC.left.val) \n right = $(BC.right.val) \n bottom = $(BC.bottom.val) \n top = $(BC.top.val) \n")
    print(" left = $(typeofBC(BC.left)) $(BC.left.val) \n right = $(typeofBC(BC.right)) $(BC.right.val) \n bottom = $(typeofBC(BC.bottom)) $(BC.bottom.val) \n top = $(typeofBC(BC.top)) $(BC.top.val) \n")
end

function typeofBC(BC)
    if is_dirichlet(BC)
        return "Dirichlet"
    elseif is_neumann(BC)
        return "Neumann"
    elseif is_navier(BC)
        return "Navier"
    elseif is_robin(BC)
        return "Robin"
    else
        return "undefined"
    end
end