% driver for multiple GRF conv plots, 2D case
% Barnett 8/15/19
clear; verb = 1;
interior = true;
curvemeth = 'n';  % see GRF_conv.m
tols = 10.^[-4:-2:-14];
figure; set(gcf,'position',[100 100 1200 900]);
for i=1:numel(tols)
  subplot(2,3,i);
  GRF_conv(tols(i),verb,curvemeth,interior);
  drawnow
end

% for paper, combine all onto a single plot (since the only thing that really
% varies is the QFS on-surf plot - label each by tol).
