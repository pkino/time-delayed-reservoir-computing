data = REC_NARMA;
data((data(:,5) == 11), 5) = NaN;
Maxeig = logspace(eigMin, eigMax, a_num);
Eiggap =  logspace(gapMax, gapMin, b_num);

horizontal = a_num;
vertical = b_num;
other = g_num;

map = zeros(vertical, horizontal,2);
for step = 1:horizontal
    for step2 = 1:vertical
        plotter = data((step-1)*vertical*other+(step2-1)*other+1:(step-1)*vertical*other+(step2-1)*other+other, 2:3);
        [map(step2,step,1), miniIndex]  = min(plotter(:,1));
        map(step2,step,2) = plotter(miniIndex,2);
    end
end

% ƒ}ƒbƒv
set(gca,'XScale','log')
X = repelem(Maxeig',vertical);
Y = repmat(Eiggap',horizontal,1);
Z = reshape(map(:,:,1),numel(map(:,:,1)),1);
c = Z;
imagesc([eigMin eigMax], [gapMax gapMin], map(:,:,1));
caxis([min(min(map(:,:,1))), 1]);
colorbar
% scatter3(X,Y,Z,700,c,'filled');
% surf(Maxeig,Eiggap,map(:,:,1));

