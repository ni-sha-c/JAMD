function [rnew,vnew]=take_half_step(r,v,F,dt,Ls)
x = r(:,1);
y = r(:,2);
z = r(:,3);

vx = v(:,1);
vy = v(:,2);
vz  = v(:,3);

ax = F(:,1);
ay = F(:,2);
az = F(:,3);

xnew = x + dt*vx + 0.5e0*dt^2*ax;
ynew = y + dt*vy + 0.5e0*dt^2*ay;
znew = z + dt*vz + 0.5e0*dt^2*az;

%Periodic BC
iselx = xnew>Ls ;
xnew(iselx) = xnew(iselx) - Ls;
isely = ynew>Ls ;
ynew(isely) = ynew(isely) - Ls;
iselz = znew>Ls ;
znew(iselz) = znew(iselz) - Ls;

iselx = xnew < 0.e0;
xnew(iselx) = xnew(iselx) + Ls;
isely = ynew < 0.e0;
ynew(isely) = ynew(isely) + Ls;
iselz = znew < 0.e0;
znew(iselz) = znew(iselz) + Ls;



rnew = [xnew, ynew, znew];


vxnew_half = vx + 0.5e0*dt*ax;
vynew_half = vy + 0.5e0*dt*ay;
vznew_half = vz + 0.5e0*dt*az;


vnew = [vxnew_half, vynew_half, vznew_half];
end