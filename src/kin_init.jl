module kin_init

meshgrid(v::AbstractVector) = meshgrid(v, v)

function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T})
    m, n = length(vy), length(vx)
    vx = reshape(vx, 1, n)
    vy = reshape(vy, m, 1)
    (repmat(vx, m, 1), repmat(vy, 1, n))
end

function meshgrid{T}(vx::AbstractVector{T}, vy::AbstractVector{T},
                     vz::AbstractVector{T})
    m, n, o = length(vy), length(vx), length(vz)
    vx = reshape(vx, 1, n, 1)
    vy = reshape(vy, m, 1, 1)
    vz = reshape(vz, 1, 1, o)
    om = ones(Int, m)
    on = ones(Int, n)
    oo = ones(Int, o)
    (vx[om, :, oo], vy[:, on, oo], vz[om, on, :])
end

#end

function initialize(Ls,Ni,Ts)

    N = Ni^3;
    DL = Ls/Ni; # initial atomic spacing 
    #Mean of MB distribution.

    ux = zeros(N,1);
    uy = zeros(N,1);
    uz = zeros(N,1);
  
    #Std of Gaussian distributions of each velocity component.
    sig = sqrt(Ts);
  
    ux = sig*rand(N,1);
    uy = sig*rand(N,1);
    uz = sig*rand(N,1);
  
    #Scaling factor for each velocity component
    scal_x = sqrt(N*Ts/sum(ux.^2));
    scal_y = sqrt(N*Ts/sum(uy.^2));
    scal_z = sqrt(N*Ts/sum(uz.^2));
  
    v = hcat(scal_x.*ux, scal_y.*uy, scal_z.*uz);
    
    #Position particles on a simple cubic lattice
    x1 = linspace(0.0 , (Ni-1)*DL, Ni);
    (x,y,z) = meshgrid(x1, x1, x1);
  
    r = hcat(x[:], y[:], z[:]);
    
    r,v
end

end


