function [Fij,Uij,rij]=force_calculation(N,r,Ls,rc2)
%Naive calculation.
    x = r(:,1);
    y = r(:,2);
    z = r(:,3);
   
   xij = spdiags(repmat(x,1,N-1),-(N-1):-1, N,N)';
   yij = spdiags(repmat(y,1,N-1),-(N-1):-1, N,N)';
   zij = spdiags(repmat(z,1,N-1),-(N-1):-1, N,N)';
   
   xj = spdiags(repmat(x,1,N-1), 1:N-1,N,N);
   yj = spdiags(repmat(y,1,N-1), 1:N-1,N,N);
   zj = spdiags(repmat(z,1,N-1), 1:N-1,N,N);
   
   xij = sparse(xij - xj);
   yij = sparse(yij - yj);
   zij = sparse(zij - zj);
   %Minimum image convention
   xij = xij - Ls*round(xij/Ls);
   yij = yij - Ls*round(yij/Ls);
   zij = zij - Ls*round(zij/Ls);
   
   rij2 = sparse(xij.^2 + yij.^2 + zij.^2);
   rij2(rij2>rc2) = 0.e0;
   [i,j,rij2] = find(rij2);
   
    %The 6 term
    rij6_inv = (1./rij2).*(1./rij2).*(1./rij2);
    %The 12 term
    rij12_inv =rij6_inv.*rij6_inv;
    
     %Force on i due to j , along rij pointing towards i.
    magfij = sparse(24.e0.*(2.e0*rij12_inv - rij6_inv)./rij2);
    
    rij6_inv = sparse(i,j,rij6_inv,N,N);
    rij12_inv = sparse(i,j,rij12_inv,N,N);
    rij2 = sparse(i,j,rij2,N,N);
    magfij = sparse(i,j,magfij,N,N);
    
    uij = sparse(4.e0.*(rij12_inv - rij6_inv));
    Uij = uij + uij';
   
    
    fijx = magfij.*xij;
    fijy = magfij.*yij;
    fijz = magfij.*zij;
    
    Fij = zeros(N,N,3);
    Fij(:,:,1) = fijx-fijx';
    Fij(:,:,2) = fijy-fijy';
    Fij(:,:,3) = fijz-fijz';
    
    rij(:,:,1) = full(xij - xij');
    rij(:,:,2) = full(yij - yij');
    rij(:,:,3) = full(zij - zij');
    
end
   

