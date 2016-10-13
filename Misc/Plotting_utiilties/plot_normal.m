figure
u1=0;
s1=0.8;
u2=0;
s2=1.6;
x=linspace(u1-3*s2,u1+3*s2,200);
y1=1/sqrt(2*pi)/s1*exp(-(x-u1).^2/(2*s1^2));
y2=1/sqrt(2*pi)/s2*exp(-(x-u2).^2/(2*s2^2));
e=0.1*log10(y2./y1);
a=[u1,u1];
b=[min(e),max(e)];
c=[u2,u2];
d=[min(e),max(e)];
subplot(2,1,1)
hold on
plot(x,y1,'-b','LineWidth',2);
plot(x,y2,':g','LineWidth',2);
plot(x,e,'.r','LineWidth',2);
axis tight
legend('Posterior density $$p(\xi|\mathbf{Y})$$','Proposal density $$\gamma(\xi|\theta)$$','$$0.1*log_{10}(\epsilon)$$',2,'Color','white') 
legend boxoff
h=0.25; v=0.2;
text('Interpreter','latex','String',['$$\epsilon = $$' num2str(10^(min(e)*10))],'Position',[((1-h)*max(x)+h*min(x)) ((1-v)*max(e)+v*min(e))],'BackgroundColor','white')
plot(a,b)
plot(c,d)

u1=-1;
s1=0.8;
u2=1;
s2=1.6;
x=linspace(u1-3*s2,u1+3*s2,200);
y1=1/sqrt(2*pi)/s1*exp(-(x-u1).^2/(2*s1^2));
y2=1/sqrt(2*pi)/s2*exp(-(x-u2).^2/(2*s2^2));
e=0.1*log10(y2./y1);
a=[u1,u1];
b=[min(e),max(e)];
c=[u2,u2];
d=[min(e),max(e)];
subplot(2,1,2)
hold on
plot(x,y1,'-b','LineWidth',2);
plot(x,y2,':g','LineWidth',2);
plot(x,e,'.r','LineWidth',2);
legend('Posterior density $$p(\xi|\mathbf{Y})$$','Proposal density $$\gamma(\xi|\theta)$$','$$0.1log_{10}(\epsilon)$$',2,'Color','white') 
legend boxoff
axis tight
h=0.25; v=0.2;
text('Interpreter','latex','String',['$$\epsilon = $$' num2str(10^(min(e)*10))],'Position',[((1-h)*max(x)+h*min(x)) ((1-v)*max(e)+v*min(e))],'BackgroundColor','white')
plot(a,b)
plot(c,d)
print_plot

figure
u1=0;
s1=1.6;
u2=0;
s2=0.8;
x=linspace(u1-3*s1,u1+3*s1,200);
y1=1/sqrt(2*pi)/s1*exp(-(x-u1).^2/(2*s1^2));
y2=1/sqrt(2*pi)/s2*exp(-(x-u2).^2/(2*s2^2));
e=0.1*log10(y2./y1);
a=[u1,u1];
b=[min(e),max(y2)];
c=[u2,u2];
d=[min(e),max(y2)];
subplot(2,1,1)
hold on
plot(x,y1,'-b','LineWidth',2);
plot(x,y2,':g','LineWidth',2);
plot(x,e,'.r','LineWidth',2);
axis tight
legend('Posterior density $$p(\xi|\mathbf{Y})$$','Proposal density $$\gamma(\xi|\theta)$$','$$0.1log_{10}(\epsilon)$$',2,'Color','white') 
legend boxoff
h=0.25; v=0.2;
text('Interpreter','latex','String',['$$\epsilon = $$' num2str(10^(min(e)*10))],'Position',[((1-h)*max(x)+h*min(x)) ((1-v)*max(e)+v*min(e))],'BackgroundColor','white')
plot(a,b)
plot(c,d)

u1=-1;
s1=1.6;
u2=1;
s2=0.8;
x=linspace(u1-3*s1,u1+3*s1,200);
y1=1/sqrt(2*pi)/s1*exp(-(x-u1).^2/(2*s1^2));
y2=1/sqrt(2*pi)/s2*exp(-(x-u2).^2/(2*s2^2));
e=0.1*log10(y2./y1);
a=[u1,u1];
b=[min(e),max(y2)];
c=[u2,u2];
d=[min(e),max(y2)];
subplot(2,1,2)
hold on
plot(x,y1,'-b','LineWidth',2);
plot(x,y2,':g','LineWidth',2);
plot(x,e,'.r','LineWidth',2);
legend('Posterior density $$p(\xi|\mathbf{Y})$$','Proposal density $$\gamma(\xi|\theta)$$','$$0.1*log_{10}(\epsilon)$$',2,'Color','white') 
legend boxoff
axis tight
h=0.25; v=0.2;
text('Interpreter','latex','String',['$$\epsilon = $$' num2str(10^(min(e)*10))],'Position',[((1-h)*max(x)+h*min(x)) ((1-v)*max(e)+v*min(e))],'BackgroundColor','white')
plot(a,b)
plot(c,d)
print_plot
