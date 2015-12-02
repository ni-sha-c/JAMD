
type Particles
    positions::Array{Float64,2}
    species::Array{Int,1}
    charges::Array{Int,1}
end

type Kinematics
    get_U::Function
    get_F::Function
end


