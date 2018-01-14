data = REC_MF_E;
data((data(:,5) == 11), 5) = NaN;
Maxeig = [-0.016, -0.008, -0.004, -0.002, -0.001,-0.0005,-0.00025, -0.000125, -0.0000625, -0.00003125];

horizontal = 10;
vertical = 20;
other = 1;
taskNum = 2;
% plotDelay = 8;
checkNum = vertical*other*taskNum;

plotter = zeros(horizontal, 3);
for step = 1:horizontal
    temp = data(plotDelay+checkNum*(step-1):taskNum*other:plotDelay+checkNum*step - taskNum, 5:6);
    [plotter(step, 1), miniIndex]  = min(temp(:,1));
    plotter(step, 2:3) = [temp(miniIndex,2), miniIndex];
end

errorbar(Maxeig, plotter(:,1), plotter(:,2),'LineWidth',2);
