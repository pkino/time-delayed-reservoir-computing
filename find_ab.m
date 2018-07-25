function [a_matrix, b_matrix] = find_ab(eigMin, eigMax, eigNum, gapMax, gapMin, gapNum)
eigen_init = eigMin; %-0.016
eigen_fin = eigMax; %3.125e-3
Max_Eigen = -logspace(eigen_init, eigen_fin, eigNum);

%% dx = -a*x + b*x_tau において、最大固有値と2番目の固有値が決まっているときにa,bを求める
tau = 80;

a_matrix = NaN(gapNum, eigNum);
b_matrix = NaN(gapNum, eigNum);
tau = 80;
dis_eig = logspace(gapMax, gapMin, gapNum); %0.005,1e-5
for step_gap = 1:gapNum
    for step_eig = 1:eigNum
        syms a b lam_sym;
        lambert = solve(lam_sym + a - b*exp(-(tau*lam_sym)) == 0, lam_sym);
        num = [0,-1];
        lambert = subs(lambert, 0 ,num);
        
        try
            [a_matrix(step_gap,step_eig), b_matrix(step_gap,step_eig)] = vpasolve([real(lambert(1))==Max_Eigen(step_eig), real(lambert(2))==Max_Eigen(step_eig)-dis_eig(step_gap)], [a, b]);
        catch
            %error('与えられた固有値配置を実現できませんでした');
        end
    end
end
end