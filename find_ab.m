function [a_matrix, b_matrix] = find_ab(eigMin, eigMax, gapMax, gapMin)
eigen_init = eigMin; %-0.016
eig_num = 4;
eigen_fin = eigMax; %3.125e-3
Max_Eigen = -logspace(eigen_init, eigen_fin, eig_num);
% for step = 1:eigen_num
%     Max_Eigen(step,1) = eigen_init*eigen_inter^(step-1);
% end

%% dx = -a*x + b*x_tau において、bと最大固有値が決まっているときにaを求める
% b_init = 0.1;
% b_max = 2;
% b_inter = 0.1;
% b_num = int32((b_max-b_init)/b_inter+1);
% b = (b_init:b_inter:b_max)';
%
% dis_eig = zeros(b_num, eigen_num);
% tau = 80;
%
% a_matrix = zeros(b_num, eigen_num);
% for step2 = 1:b_num
%     for step3 = 1:eigen_num
%         syms a
%         a_matrix(step2,step3) = double(solve(Max_Eigen(step3) + a -b(step2)*exp(-(tau*Max_Eigen(step3))) == 0, a));
%
%         syms lam_sym;
%         a_lam = a_matrix(step2,step3);
%         b_lam = b(step2);
%         lambert = solve(lam_sym + a_lam -b_lam*exp(-(tau*lam_sym)) == 0, lam_sym);
%         num = [0,1];
%         lambert = subs(lambert, 0 ,num);
%         dis_eig(step2,step3) = lambert(1) - real(lambert(2));
%     end
% end

%% dx = -a*x + b*x_tau において、最大固有値と2番目の固有値が決まっているときにa,bを求める
gap_num = 4;
tau = 80;

a_matrix = NaN(gap_num, eig_num);
b_matrix = NaN(gap_num, eig_num);
tau = 80;
dis_eig = logspace(gapMax, gapMin, gap_num); %0.005,1e-5
for step_gap = 1:gap_num
    for step_eig = 1:eig_num
        syms a b lam_sym;
        lambert = solve(lam_sym + a -b*exp(-(tau*lam_sym)) == 0, lam_sym);
        num = [0,-1];
        lambert = subs(lambert, 0 ,num);
        
        try
            [a_matrix(step_gap,step_eig), b_matrix(step_gap,step_eig)] = vpasolve([real(lambert(1))==Max_Eigen(step_eig), real(lambert(2))==Max_Eigen(step_eig)-dis_eig(step_gap)], [a, b]);
        catch
        end
    end
end
end