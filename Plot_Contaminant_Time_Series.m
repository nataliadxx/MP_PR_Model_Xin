%----------- Plot normalized contaminant series
figure(3)
clf
subplot (3,1,1)
plot (tt(1:Nm), x(1:Nm)/x(1), 'b-')
ylabel ('x/x_o','fontweight','bold','fontsize',10)

subplot (3,1,2)
plot (tt(1:Nm),UPx(1:Nm),'g-')
hold on
ylabel ('UP(gm/m^2/d)','fontweight','bold','fontsize',10)

subplot (3,1,3)
plot (tt(1:Nm),LEx(1:Nm),'r--')
ylabel ('LE (gm/m^2/d)','fontweight','bold','fontsize',10)
xlabel ('Time (d)','fontweight','bold','fontsize',10)