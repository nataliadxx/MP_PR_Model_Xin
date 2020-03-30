%----------- Plot normalized contaminant series
figure(3)
clf
subplot (2,1,1)
plot (tt(1:Nm), x(1:Nm)/x(1), 'b-')
ylabel ('x/x_o','fontweight','normal','fontsize',15)

subplot (2,1,2)
plot (tt(1:Nm),UPx(1:Nm),'g-','linewidth',2)
hold on
plot (tt(1:Nm),LEx(1:Nm),'r--','linewidth',2)
legend('UP','LE')
xlabel ('Time (d)','fontweight','normal','fontsize',15)
ylabel ('UP & LE (g mm^{-2} d^{-1})','fontweight','normal','fontsize',15)

