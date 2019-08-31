% driver for multiple GRF conv plots, 2D case, int/ext Laplace
% Barnett 8/15/19
clear; verb = 1;
interior = true;
curvemeth = 'n';                % 'i' or 'n'; see GRF_conv.m
tols = 10.^[-4:-2:-14];
figure; set(gcf,'position',[100 100 1200 900]);
for i=1:numel(tols)
  subplot(2,3,i);
  [Ns, es{i}, cns{i}] = GRF_conv(tols(i),verb,curvemeth,interior);
  es{i} = abs(es{i});
  drawnow
end

% for paper, combine all onto a single plot (since the only thing that really
% varies is the QFS on-surf plot - label each by tol).
figure;
hn = semilogy(Ns,max(es{1}(4,:),eps),'b.-');   % native far
hold on; hk = plot(Ns,es{1}(2,:),'go-');   % Kress
for i=1:numel(tols)
  hq = plot(Ns, es{i}(1,:),'k+-');     % QFS on-surf
  hf = plot(Ns, es{i}(3,:),'rd-');     % QFS far
  text(Ns(end-1)-30,3*es{i}(1,end-1),sprintf('$\\epsilon=10^{%d}$',log10(tols(i))),'interpreter','latex');
end
legend([hn hk hf hq],'native (far)','Kress (on-surf. L^2)', 'QFS (far)','QFS (on-surf. L^2)','location','north');
title(sprintf('GRF conv, int=%d, curvemeth=%s',interior,curvemeth));  % *** edit
xlabel('N'); ylabel('error');
axis([min(Ns),max(Ns),eps,1]);

