function [r,v]=initialize(Ls,Ni,Ts)

N = Ni^3;
DL = Ls/Ni; % initial atomic spacing 
%Mean of MB distribution.
  ux = zeros(N,1);
  uy = zeros(N,1);
  uz = zeros(N,1);
  
  %Std of Gaussian distributions of each velocity component.
  sig = repmat(sqrt(Ts), N, 1);
  
  ux = normrnd(ux,sig);
  uy = normrnd(uy,sig);
  uz = normrnd(uz,sig);
  %Scaling factor for each velocity component
  scal_x = sqrt(N*Ts/sum(ux.^2));
  scal_y = sqrt(N*Ts/sum(uy.^2));
  scal_z = sqrt(N*Ts/sum(uz.^2));
  
  v = [scal_x.*ux, scal_y.*uy, scal_z.*uz];
  %Position particles on a simple cubic lattice
  x1 = linspace(0.0 , (Ni-1)*DL, Ni);
  [x,y,z] = meshgrid(x1, x1, x1);
  
  r = [x(:), y(:), z(:)];
end