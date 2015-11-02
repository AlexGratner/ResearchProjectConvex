%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static optimization for planetary gear                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cvx_begin gp
    variables n rg1 bg1 n1 rg2 bg2 n2 
    minimize(pi * power(rg1,2) * 1.5 * bg1 + pi * power(rg2,2) *1.5 *bg2)
    subject to       
        %Planetary gear
        n1 * power(rg1,-2) * power(bg1,-1) + (power(n1,-1)+ 2*power(n1,-2)+4*power(n1,-3)+8*power(n1,-4)+16*power(n1,-5)+32*power(n1,-6)+64*power(n1,-7))*power(rg1,-2)*power(bg1,-1) <= k_inv;
        n2 * power(rg2,-2) * power(bg2,-1)*power(n1,-1) + (power(n2,-1)+ 2*power(n2,-2)+4*power(n2,-3)+8*power(n2,-4)+16*power(n2,-5)+32*power(n2,-6)+64*power(n2,-7))*power(rg2,-2)*power(bg2,-1)*power(n1,-1) <= k_inv;
        rg1 == rg2;
        1/20 <= bg1/rg1 <= 4;
        1/20 <= bg2/rg2 <= 4;
        3 <= n1;
        3 <= n2;
        n1*n2*power(n,-1) ==1;
        n == 60;
cvx_end