%----------- Plot degree of saturation time series
figure(2)
clf
subplot (4,1,1)
plot (tt(1:Nm), Is, 'b-')
ylabel ('P (mm/d)','fontweight','bold','fontsize',10)

subplot (4,1,2)
plot (tt(1:Nm),s(1:Nm),'k-')
ylabel ('s','fontweight','bold','fontsize',10)

subplot (4,1,3)
plot (tt(1:Nm),ET(1:Nm),'g:','linewidth',3)
ylabel ('ET (mm/d)','fontweight','bold','fontsize',10)

subplot (4,1,4)
plot (tt(1:Nm),LQ(1:Nm),'r--')
xlabel ('Time (d)','fontweight','bold','fontsize',10)
ylabel ('LQ (mm/d)','fontweight','bold','fontsize',10)

