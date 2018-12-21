data = NRMSE;
Maxeig = logspace(eigMin, eigMax, eigNum);
Eiggap =  logspace(gapMax, gapMin, gapNum);

% set indexes of saved data
NRMSEIndex = 6;
paramsStartIndex = 1; % the start row of the saved parameters
paramsLastIndex = 5;
cIndex = 3;
if NRMSEIndex > paramsLastIndex
    paramsLastIndex = NRMSEIndex + 2;
end

[bestPerform, bestIndex] = min(data(:,1));
bestC = data(bestIndex, cIndex);
bestParams = data(bestIndex, paramsStartIndex:paramsLastIndex)

horizontal = eigNum;
vertical = gapNum;
other = gammaNum*cNum;

map = zeros(vertical, horizontal,2);
for step = 1:vertical
    for step2 = 1:horizontal
        plotter = data((step-1)*vertical*other+(step2-1)*other+1:(step-1)*vertical*other+(step2-1)*other+other, 1:paramsLastIndex);
        
        % fixed c value
        plotter = plotter(plotter(:,cIndex) == bestC,:);
        
        
        [map(step,step2,1), miniIndex]  = min(plotter(:,NRMSEIndex));
        %         miniIndex
        map(step,step2,2) = plotter(miniIndex,end);
        
        % If you want to know each value, uncomment.
        %         miniIndex
        %         c=plotter(miniIndex,cIndex)
    end
end

%% draw a map
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


