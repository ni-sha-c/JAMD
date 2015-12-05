
module MD_calc

#################################

function take_a_step(N,r,v,F,t,dt,Ls,STEPS_thermostat,alpha,Ts)

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

function force_calculation(N,r,Ls,rc2)

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
    Cut = rij2_mat.<=rc2;
    
    rij6i_mat = (1./rij2_mat.^3);
    rij12i_mat = (rij6i_mat.^2);
    
    Ftemp = (24./rij2_mat.*(2.*rij12i_mat-rij6i_mat));
    
    for i = 1:N
        Ftemp[i,i] = 0;
        rij6i_mat[i,i] = 0;
        rij12i_mat[i,i] = 0;
    end

    Fij=zeros(N,N,3)

    Fij[:,:,1] = Ftemp.*rij_mat[:,:,1].*Cut;
    Fij[:,:,2] = Ftemp.*rij_mat[:,:,2].*Cut;
    Fij[:,:,3] = Ftemp.*rij_mat[:,:,3].*Cut;
    Uij = 4*(rij12i_mat-rij6i_mat).*Cut;

    rij = rij_mat;
    
    Fij,Uij,rij

end

end
