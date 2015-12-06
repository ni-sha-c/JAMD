function [P,U]=P_and_U(N,Vs,rij,Fij,Uij,T_inst,P_LRC,E_LRC)    

    fijx = triu(Fij(:,:,1));
    fijy = triu(Fij(:,:,2));
    fijz = triu(Fij(:,:,3));

    xij = triu(rij(:,:,1));
    yij = triu(rij(:,:,2));
    zij = triu(rij(:,:,3));

    rdotf =  fijx.*xij + fijy.*yij + fijz.*zij;
    rdotf = rdotf(:);
    w = 1.e0/3.e0*sum(rdotf);

    pterm1 = w/Vs;
    pterm2 = N/Vs*T_inst;
    
    P = pterm1 + pterm2;
    uij = triu(Uij);
    uij = uij(:);
    U = sum(uij)/N;
    
    P =P + P_LRC;
    U =U+E_LRC;

end
