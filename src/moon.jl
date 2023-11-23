function calc_proj_dist2(p1, p2)
    """
    determines distance between moon and cell on solar grid 

    p1: vector from observer to cell on solar grid
    p2: vector from observer to moon 
    """
    x1 = p1[1]
    x2 = p2[1] 
    y1 = atan(p1[2] / x1) 
    y2 = atan(p2[2] / x2)
    z1 = atan(p1[3] / x1)
    z2 = atan(p2[3] / x2)
    return (y1 - y2)^2 + (z1 - z2)^2
end

function quad_limb_darkening_optical(μ::T) where T
    """
    limb darkening prescription for optical based on mu angle  
    """
    μ < zero(T) && return 0.0
    return 0.28392 + 1.36896*μ - 1.75998*μ^2 + 2.22154*μ^3 - 1.56074*μ^4 + 0.44630*μ^5
end

function quad_limb_darkening_NIR(μ::T) where T
    """
    limb darkening prescription for NIR based on mu angle  
    """
    μ < zero(T) && return 0.0
    return 0.59045 + 1.41938*μ - 3.01866*μ^2 + 3.99843*μ^3 - 2.67727*μ^4 + 0.068758*μ^5
end