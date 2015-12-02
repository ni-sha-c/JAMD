module particle_manipulations



global num_atoms = Int64
global positions = Array{Float64,2}
global species = ""

function read(filename)
    global num_atoms, species, positions
	f = open(filename, "r") #Read a standard .xyz file type

    i = 1
	for l in enumerate(eachline(f))

       data = readdlm(IOBuffer(l[2]))
       
	   if(l[1]==1)
			num_atoms = convert(Int64,data[1])
	   		positions = zeros(num_atoms,3)
			 
		end
	   if(l[1]>2)
	
	   species = string(species," ", convert(ASCIIString,data[1])) 
	   positions[i,1] = convert(Float64,data[2]) 
	   positions[i,2] = convert(Float64,data[3]) 
	   positions[i,3] = convert(Float64,data[4])
	   i = i + 1
	   end
    end

    close(f) 
end

function check_same_particle(p1::Particles, p2::Particles)
    if length(p1)!=length(p2)
        return false
    elseif any(p1.positions!=p2.positions)
        return false
    elseif any(p1.species!=p2.positions)
        return false
    elseif any(p1.charges!=p1.charges)
        return false
    else
        return true
    end
end

function get_R(p::Particles, i::Int, j::Int)
    #dist  = p.positions[i,1:] - p.positions[j,1:]
    #sqrt(sum(dist.^2))
end

end
