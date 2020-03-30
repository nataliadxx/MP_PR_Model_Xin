%----------- Plot degree of saturation time series
figure(2)
clf

subplot (2,1,1)
plot (tt(1:Nm),s(1:Nm),'k-')
ylabel ('s','fontweight','normal','fontsize',15)

subplot (2,1,2)
plot (tt(1:Nm),ET(1:Nm),'g:','linewidth',2)
hold on
plot (tt(1:Nm),LQ(1:Nm),'r--','linewidth',2)
legend('ET','LQ')
xlabel ('Time (d)','fontweight','normal','fontsize',15)
ylabel ('ET & LQ (mm d^{-1})','fontweight','normal','fontsize',15)