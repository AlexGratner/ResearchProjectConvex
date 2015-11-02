%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static optimization for motor                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


cvx_begin gp
    variables l r n 
    minimize( pi * power(r,2) * l)
    subject to
        %Motor
        inertia_ratio^2/Jload^2*Cmj^2*Trms^2 * power(r,3)*power(n,2) + 2*inertia_ratio/Jload*Cmj*(Trms^2) * power(r,-1)*power(l,-1)+Trms^2*power(n,-2)*power(l,-2)*power(r,-5) <= Cm^2;
        l <= 0.18;
        9 <= n <= 100;
        0.03 <= r <= 0.07;
        lb <= l/r <= ub;

cvx_end