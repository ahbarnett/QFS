% driver for multiple GRF conv plots
% Barnett 8/15/19
clear; verb = 0;
tols = 10.^[-4:-2:-14];
figure(1); clf; set(gcf,'position',[100 100 1200 900]);
for i=1:numel(tols)
  subplot(2,3,i);
  GRF_conv(tols(i),verb);
  drawnow
end
