
module MD_calc

function take_half_step(r,v,F,dt,Ls)
    
    x=r[:,1];
    y=r[:,2];
    z=r[:,3];
    
    vx=v[:,1];
    vy=v[:,2];
    vz=v[:,3];
    
    ax=F[:,1];
    ay=F[:,2];
    az=F[:,3];
    
    xnew = x + dt*vx + 0.5*dt^2*ax;
    ynew = y + dt*vy + 0.5*dt^2*ay;
    znew = z + dt*vz + 0.5*dt^2*az;

    #Periodic BC   
    xnew-=Ls.*(xnew.>Ls)
    ynew-=Ls.*(ynew.>Ls)
    znew-=Ls.*(znew.>Ls)

    xnew+=Ls.*(xnew.<0)
    ynew+=Ls.*(ynew.<0)
    znew+=Ls.*(znew.<0)

    rnew = hcat(xnew, ynew, znew);

    vxnew_half = vx + 0.5*dt*ax;
    vynew_half = vy + 0.5*dt*ay;
    vznew_half = vz + 0.5*dt*az;

    vnew = hcat(vxnew_half, vynew_half, vznew_half);

    rnew,vnew
end

function take_one_step(N,v,F,t,dt,STEPS_thermostat,alpha,Ts)
   
    vx = v[:,1];
    vy = v[:,2];
    vz  = v[:,3];

    ax = F[:,1];
    ay = F[:,2];
    az = F[:,3];

    vxnew = vx + 0.5*dt*ax;
    vynew = vy + 0.5*dt*ay;
    vznew = vz + 0.5*dt*az;

    vnew = hcat(vxnew, vynew, vznew);
    
    c = [sum(vnew[:,1]) sum(vnew[:,2]) sum(vnew[:,3])]/N;
    v_for_temp = broadcast(-,vnew,c);
 
    #Calculate Instantaneous Temperature
    T_inst = 1./(3.*N)*sum(v_for_temp.^2);

    #Velocity Rescaling
    if(t<STEPS_thermostat) 
        v_scal = sqrt(1 + alpha*(Ts/T_inst - 1.));
        vnew = v_scal*vnew;
    end
    
    vnew,T_inst
end

#################################

function take_one_step_good(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts)

vnew = v+dt*F;
v_c = 0.5*(vnew+v);

c = [sum(v_c[:,1]) sum(v_c[:,2]) sum(v_c[:,3])]/N;
SUM = sum((v_c[:,1]-c[1]).^2+(v_c[:,2]-c[2]).^2+(v_c[:,3]-c[3]).^2);
T_inst = 2/(3*N)*0.5*SUM;

if t < STEPS_thermostat
vnew = vnew*sqrt( 1+alpha*(Ts/T_inst-1) );
end

rnew = r+dt*vnew;

Pos2 = rnew.>=Ls;
Neg2 = rnew.<0;
rnew = rnew-Ls.*Pos2+Ls.*Neg2;

rnew,vnew,T_inst

end

#################################

function force_calculation_mod(N,r,Ls,rc2)

    r_matrix_x = repmat(r[:,1],1,N)-repmat(r[:,1],1,N)';
    r_matrix_y = repmat(r[:,2],1,N)-repmat(r[:,2],1,N)';
    r_matrix_z = repmat(r[:,3],1,N)-repmat(r[:,3],1,N)';
    Pos = r_matrix_x.>Ls/2.;
    Neg = r_matrix_x.<-Ls/2.;
    r_matrix_x = r_matrix_x-Ls.*(Pos-Neg);
    Pos = r_matrix_y.>Ls/2.;
    Neg = r_matrix_y.<-Ls/2.;
    r_matrix_y = r_matrix_y-Ls.*(Pos-Neg);
    Pos = r_matrix_z.>Ls/2.;
    Neg = r_matrix_z.<-Ls/2.;
    r_matrix_z = r_matrix_z-Ls.*(Pos-Neg);

    rij_mat=zeros(N,N,3)

    rij_mat[:,:,1] = r_matrix_x;
    rij_mat[:,:,2] = r_matrix_y;
    rij_mat[:,:,3] = r_matrix_z;     
        
    rij2_mat = r_matrix_x.*r_matrix_x+r_matrix_y.*r_matrix_y+r_matrix_z.*r_matrix_z;
#    Cut = [rij2_mat<=rc2];
    
    rij6i_mat = (1./rij2_mat.^3);
    rij12i_mat = (rij6i_mat.^2);
    
    Ftemp = (24./rij2_mat.*(2.*rij12i_mat-rij6i_mat));
    
    for i = 1:N
        Ftemp[i,i] = 0;
        rij6i_mat[i,i] = 0;
        rij12i_mat[i,i] = 0;
    end

    Fij=zeros(N,N,3)

    Fij[:,:,1] = Ftemp.*rij_mat[:,:,1]#.*Cut;
    Fij[:,:,2] = Ftemp.*rij_mat[:,:,2]#.*Cut;
    Fij[:,:,3] = Ftemp.*rij_mat[:,:,3]#.*Cut;
    Uij = 4*(rij12i_mat-rij6i_mat)#.*Cut;

    rij = rij_mat;
    
    Fij,Uij,rij

end

function force(N,r,Ls,rc2)

x = r[:,1];
y = r[:,2];
z = r[:,3];

xij = diag(repmat(x,1,N-1),-(N-1):-1, N,N)';
yij = diag(repmat(y,1,N-1),-(N-1):-1, N,N)';
zij = diag(repmat(z,1,N-1),-(N-1):-1, N,N)';

xj = diag(repmat(x,1,N-1), 1:N-1,N,N);
yj = diag(repmat(y,1,N-1), 1:N-1,N,N);
zj = diag(repmat(z,1,N-1), 1:N-1,N,N);

xij = sparse(xij - xj);
yij = sparse(yij - yj);
zij = sparse(zij - zj);

#Minimum image convention
xij = xij - Ls*round(xij/Ls);
yij = yij - Ls*round(yij/Ls);
zij = zij - Ls*round(zij/Ls);

#rij2 = sparse(xij.^2 + yij.^2 + zij.^2);
#rij2(rij2>rc2) = 0.e0; # TOOK OUT CUTOFF FOR NOW
(i,j,rij2) = find(rij2);

#The 6 term
rij6_inv = (1./rij2).*(1./rij2).*(1./rij2);
#The 12 term
rij12_inv =rij6_inv.*rij6_inv;

#Force on i due to j, along rij pointing towards i.
magfij = sparse(24.e0.*(2.e0*rij12_inv - rij6_inv)./rij2);

rij6_inv = sparse(i,j,rij6_inv,N,N);
rij12_inv = sparse(i,j,rij12_inv,N,N);
rij2 = sparse(i,j,rij2,N,N);
magfij = sparse(i,j,magfij,N,N);

uij = sparse(4.*(rij12_inv - rij6_inv));
Uij = uij + uij';


fijx = magfij.*xij;
fijy = magfij.*yij;
fijz = magfij.*zij;

Fij = zeros(N,N,3);
Fij[:,:,1] = fijx-fijx';
Fij[:,:,2] = fijy-fijy';
Fij[:,:,3] = fijz-fijz';

rij[:,:,1] = full(xij - xij');
rij[:,:,2] = full(yij - yij');
rij[:,:,3] = full(zij - zij');

Fij,Uij,rij

end


end
