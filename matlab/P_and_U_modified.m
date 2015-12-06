function [P,u]=P_and_U_modified(N,Vs,r,F,U,T_inst,P_LRC,E_LRC)    
    u = sum(U)/N;
    P = 0.e0;
   % P =P + P_LRC;
    u =u+E_LRC;

end
