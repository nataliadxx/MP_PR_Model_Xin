%------------Plot LAI, photosynthesis and root depth
figure(4)
clf
subplot(2,1,1)
plot(tt(1:Nm),LAI(1:Nm),'g-','linewidth',2)
ylabel('LAI','fontweight','normal','fontsize',15)

subplot(2,1,2)
plot(tt(1:Nm),Zr(1:Nm),'b-')
ylabel('Root depth (mm)','fontweight','normal','fontsize',15)
xlabel('Time (d)','fontweight','normal','fontsize',15)