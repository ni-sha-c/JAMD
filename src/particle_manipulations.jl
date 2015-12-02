
function read(filename)
    f = open(filename, "r") #Read a standard .xyz file type

    num_atoms = int(readline(f))
    throwaway = chomp(readline(f)) 
    
    positions = zeros(Float64, num_atoms, 3)
    species = zeros(Int, num_atoms)
    charges = zeros(Int, num_atoms)


    for (line,i) in enumerate(each_line(f))
        fields = split(line)
        positions[i,:] = map(float, fields[2:])
    end

    close(f)

    #Make Particles 
end

function ==(p1::Particles, p2::Particles)
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
    dist  = p.positions[i,1:] - p.positions[j,1:]
    sqrt(sum(dist.^2))
end
