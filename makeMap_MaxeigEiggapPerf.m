data = NRMSE;
Maxeig = logspace(eigMin, eigMax, eigNum);
Eiggap =  logspace(gapMax, gapMin, gapNum);
saveNum = 5;

cIndex = 3;
[bestPerform, bestIndex] = min(data(:,saveNum+1));
bestC = data(bestIndex, cIndex);
bestParams = data(bestIndex, 1:saveNum)

horizontal = eigNum;
vertical = gapNum;
other = gammaNum*cNum;

map = zeros(vertical, horizontal,2);
for step = 1:vertical
    for step2 = 1:horizontal
        plotter = data((step-1)*vertical*other+(step2-1)*other+1:(step-1)*vertical*other+(step2-1)*other+other, 1:saveNum+2);
        
        % c固定
        plotter = plotter(plotter(:,cIndex) == bestC,:);
        
        
        [map(step,step2,1), miniIndex]  = min(plotter(:,end-1));
        map(step,step2,2) = plotter(miniIndex,end);
        
        % 各値が知りたいときはコメントアウトを外す
        %         miniIndex
        %         c=plotter(miniIndex,3)
    end
end

%% マップ描画
figure;

set(gca,'XScale','log')
X = repmat(Eiggap',horizontal,1);
Y = repelem(Maxeig',vertical);
Z = reshape(map(:,:,1),numel(map(:,:,1)),1);
c = Z;
imagesc([eigMin eigMax], [gapMax gapMin], map(:,:,1));
caxis([min(min(map(:,:,1))), 1]);
colorbar;

xlabel('maximum eigenvalue'); ylabel('gap');

az = 180;
el = 90;
view(az, el);


