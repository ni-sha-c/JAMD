function [vnew,T_inst]=take_one_step(N,v,F,t,dt,STEPS_thermostat,alpha,Ts)
vx = v(:,1);
vy = v(:,2);
vz  = v(:,3);

ax = F(:,1);
ay = F(:,2);
az = F(:,3);

vxnew = vx + 0.5e0*dt*ax;
vynew = vy + 0.5e0*dt*ay;
vznew = vz + 0.5e0*dt*az;


vnew = [vxnew, vynew, vznew];
 c = sum(vnew)/N;
 c = repmat(c,N,1);
 %Calculate Instantaneous Temperature
T_inst = 1.e0/3.e0/N*sum(sum((vnew-c).^2));


%Velocity Rescaling.
if(t<STEPS_thermostat) 
 v_scal = sqrt(1 + alpha*(Ts/T_inst - 1.e0));
  vnew = repmat(v_scal,N,3).*vnew;

end
end