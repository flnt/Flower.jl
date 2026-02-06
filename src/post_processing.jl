"""
TODO constant mesh spacing
"""
function calculate_centroid(x, y, volume_cell)
    # x and y are scalar node coord

    # println(y)

    # println(volume_cell)

    integral_x = sum(x .* volume_cell)
    integral_y = sum(y .* volume_cell)
    integral_1 = sum(volume_cell)
    if integral_1 > 0.0
        x_c = integral_x / integral_1
        y_c = integral_y / integral_1
    else
        x_c = 0.0
        y_c = 0.0
        println("\033[31m error calculate_centroid\033[0m")
    end

    return (x_c, y_c)
end

function calculate_circularity(perimeter_bubble, area)
    # area = pi r 2
    # r = sqrt(area/pi)
    # perim= 2 sqrt(area*pi)
    # perim = 2pi r
    # area = volume_fraction * dx dy ou dcap
    if perimeter_bubble > 0.0
        perimeter_circle = 2 * sqrt(pi * area)
        circularity = perimeter_circle / perimeter_bubble
    else
        circularity = 0.0
        println("\033[31m error calculate_circularity\033[0m")
    end
    return circularity
end

function calculate_rise_velocity(v, volume_cell)
    # v velocity
    integral_u = sum(v .* volume_cell)
    integral_1 = sum(volume_cell)
    if integral_1 > 0.0
        U_c = integral_u / integral_1
    else
        U_c = 0.0
        println("\033[31m error calculate_rise_velocity\033[0m")
    end
    return U_c
end