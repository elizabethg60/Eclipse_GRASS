function calc_proj_dist2(p1, p2)
    """
    determines distance in radians between moon and cell on solar grid

    p1: vector from observer to cell on solar grid
    p2: vector from observer to moon 
    """
    return acos(dot(p1, p2)/(norm(p1)*norm(p2)))
end

function quad_limb_darkening_optical(μ::T, index, zenith_angle) where T
    """
    limb darkening prescription for optical based on mu angle  
    """
    μ < zero(T) && return 0.0
    I = 0.28392 + 1.36896*μ - 1.75998*μ^2 + 2.22154*μ^3 - 1.56074*μ^4 + 0.44630*μ^5
    cell_airmass = 1/cosd(zenith_angle)
    ext_factor =  -(cell_airmass*ext_coef[index])/2.5 
    return I # * 10^ext_factor #being calculated separated in compute_rv (will have to comment out for other locations)
end

function quad_limb_darkening_NIR(μ::T) where T
    """
    limb darkening prescription for NIR based on mu angle  
    """
    μ < zero(T) && return 0.0
    return 0.59045 + 1.41938*μ - 3.01866*μ^2 + 3.99843*μ^3 - 2.67727*μ^4 + 0.068758*μ^5
end
