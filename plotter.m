% Cross-track
figure
xx = ceil(length(x)/2);
a1 = plot(y*lrad*1e9,ones(size(h_lube_disk(xx,:,1)))*h0*1e9,'b',y*lrad*1e9,fh*h0*1e9-0.3*ones(size(h_lube_slider(xx,:,1)))*h0*1e9,'b',...
    y*lrad*1e9,h_lube_disk(xx,:,round(tf/dt/s_interval*1/3))*h0*1e9,'-.r',y*lrad*1e9,fh*h0*1e9-h_lube_slider(xx,:,round(tf/dt/s_interval*1/3))*h0*1e9,'-.r',...
    y*lrad*1e9,h_lube_disk(xx,:,round(tf/dt/s_interval*2/3))*h0*1e9,'-^k',y*lrad*1e9,fh*h0*1e9-h_lube_slider(xx,:,round(tf/dt/s_interval*2/3))*h0*1e9,'-^k',...
    y*lrad*1e9,h_lube_disk(xx,:,round(tf/dt/s_interval))*h0*1e9,':G',y*lrad*1e9,fh*h0*1e9-h_lube_slider(xx,:,round(tf/dt/s_interval))*h0*1e9,':G','markers',1.5);
xlim([-40,40])
ylim([0, 4])
xticks([-40:10:40])
yticks([0:0.5:4])
xlabel('Cross-track (Y) direction (nm)')
ylabel('Normal Co-ordinate (Z) (nm)')
set(a1,'linewidth',2);
grid on
lgd = legend([a1(1) a1(3) a1(5) a1(7)],'time = 0','time = 0.33tf', 'time = 0.67tf', 'time = tf');
set(gca,'fontsize', 14)

% Downtrack
figure
yy = ceil(length(y)/2);
a1 = plot(x*lrad*1e6,h_lube_disk(:,yy,end)*h0*1e9,':g',x*lrad*1e6,fh*h0*1e9-h_lube_slider(:,yy,end)*h0*1e9,':g','markers',1.5);
ylim([0, 4])
yticks([0:0.5:4])
xlabel('Down-track (X) direction (\mum)')
ylabel('Normal Co-ordinate (Z) (nm)')
set(a1,'linewidth',1.5);
grid on
lgd = legend('time = tf');
set(gca,'fontsize', 14)