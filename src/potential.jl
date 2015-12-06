module potential 
include("Definitions.jl")
function LJ(sigma::Float64,epsilon::Float64,p::Particles)
   
	rc = 2.5
    N = length(p)
    E = 0.0
    F = zeros(Float64, N, 3)

    for i = 1:N
        for j = i+1:N
            r = get_R(p, i, j)
            E += 4*epsilon*((sigma/r)^12.0 - (sigma/r)^6.0)
            F = -24*epsilon*(2*(sigma^12.0/r^13.0) - (sigma^6.0/r^7.0))

            diff = p.positions[i,:] - p.positions[j,:] 
            F[i,:] -= f*diff/r
            F[j,:] += f*diff/r
        end
    end
end

end
