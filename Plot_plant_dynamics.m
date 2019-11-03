%------------Plot LAI, photosynthesis and root depth
figure(4)
clf
subplot(3,1,1)
plot(tt(1:Nm),LAI(1:Nm),'g-')
ylabel('LAI','fontweight','normal','fontsize',10)

subplot(3,1,2)
plot(tt(1:Nm),An(1:Nm),'r-')
ylabel('Photosynthesis (kg C/d)','fontweight','normal','fontsize',10)

subplot(3,1,3)
plot(tt(1:Nm),Zr(1:Nm),'b-')
ylabel('Root depth (mm)','fontweight','normal','fontsize',10)