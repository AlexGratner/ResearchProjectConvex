in_for_pole = 'What is the pole values (expecting negative values)? ';
poles = -input(in_for_pole);

lb_k = damp*min(poles);                                                     %Lower bound for k
ub_delay = 1/poles(2) + 1/poles(2) + 1/poles(3) + 1/poles(4);               %Upper bound for delay time             
lb_delay = ub_delay - damp/lb_k;                                            %Lower bound for delay time

in_for_delay = ['What is the reaction delay time expected (should be in between of ',num2str(lb_delay), ' and ' num2str(ub_delay), ')? '];
disp(in_for_delay);
t_delay = input('Expecting positive');

cvx_begin gp
    variables l r n rg1 bg1 n1 rg2 bg2 n2 kl rl ll
    minimize( pi * power(r,2) * l + pi * power(rg1,2) * 1.5 * bg1 + pi * power(rg2,2) *1.5 *bg2 + pi * power(rl,2) * ll)
    subject to
        %Motor
        inertia_ratio^2/Jload^2*Cmj^2*Trms^2 * power(r,3)*power(n,2) + 2*inertia_ratio/Jload*Cmj*(Trms^2) * power(r,-1)*power(l,-1)+Trms^2*power(n,-2)*power(l,-2)*power(r,-5) <= Cm^2;
        l <= 0.18;
        9 <= n <= 100;
        0.03 <= r <= 0.07;
        lb <= l/r <= ub;
        
        %Planetary gear
        n1 * power(rg1,-2) * power(bg1,-1) + (power(n1,-1)+ 2*power(n1,-2)+4*power(n1,-3)+8*power(n1,-4)+16*power(n1,-5)+32*power(n1,-6)+64*power(n1,-7))*power(rg1,-2)*power(bg1,-1) <= k_inv;
        n2 * power(rg2,-2) * power(bg2,-1)*power(n1,-1) + (power(n2,-1)+ 2*power(n2,-2)+4*power(n2,-3)+8*power(n2,-4)+16*power(n2,-5)+32*power(n2,-6)+64*power(n2,-7))*power(rg2,-2)*power(bg2,-1)*power(n1,-1) <= k_inv;
        rg1 == rg2;
        1/20 <= bg1/rg1 <= 4;
        1/20 <= bg2/rg2 <= 4;
        3 <= n1;
        3 <= n2;
%         0.03 <= rg1 <= 0.1;
        n1*n2*power(n,-1) ==1;
        
        % Shaft
        lb_k <= kl;
        ub_delay - t_delay <= damp*power(kl,-1);
        kl == G*pi/2*power(rl,4)*power(ll,-1);
        0.05 <= ll;
        
        
        
        
%         kl <= 1.88e4;
%         kl*ll*power(rl,-4)*2/pi/G == 1;
%         0.03 <= ll <= 0.06;
%         0.003<= rl;
cvx_end