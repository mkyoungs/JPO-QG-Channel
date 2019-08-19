% clear
% clear global
global nx ny f0 r kappa
scale = 2;
nx=128*scale;ny=64*scale;


Cr('ch2.h','init.c');
Cr('#define','NX',nx);
Cr('#define','NY',ny);
Cr('#define','NYP',ny+1);
Cr('#define','SHIFT',round(log2(nx)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


r1=0;Cr('float','r1',r1);
r=0.03;Cr('float','r',r); % 0.03 in MITgcm
%r1=5e-3;Cr('float','r1',r1); % same piston velocity as MITgcm 5e-4

H0 = -2; % overturning in Sv
tau0m = 0.1; % in N/m^2 kg/m/s^2
tau0c = tau0m;
top0 = 1;% topography in km

dtinv=64;Cr('float','dtinv',dtinv);
tmax=120000/dtinv;Cr('float','tmax',tmax);
kappa=1.3;Cr('float','kappa',kappa);%1.3
del=0.2;Cr('float','del',del); %0.2
rd=15;Cr('float','rd',rd);%15

% this leads to a density difference of 0.4 kg/m^3
% f0^2/fr1/H1*1000/(3600*3600*24*24)/10*1000
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
frm=0;



width=2000;Cr('float','width',width);
%rd=38.8;Cr('float','rd',rd);


f0=-8.64;Cr('float','f0',f0); % -8.64
beta=1*8.64e-4;Cr('float','beta',beta); % 8.64e-4 1/day/km

totaldepth = 4;
H2 = totaldepth/(1+del);
H1 = totaldepth - H2;
fr1=1.0/rd/rd/(1.0+del);fr2=del*fr1; % units of 1/km/km



%kappa=8.64;Cr('float','kappa',kappa);



tpl=dtinv*160/dtinv;Cr('int','tpl',tpl);

dx=width/ny;
nyp=ny+1;

x=[0.5:1:nx]*dx;
xp=[0:nx]*dx;
y=[0.5:1:ny]'*dx-width/2;
yp=[0:ny]'*dx-width/2;
yy=yp*ones(1,nx);




%Cr('float','Hy',Hy);

%tau=tau0*1000*3600^2*24^2*cos(pi*yp/width); % in km/day units
taumean = tau0m*2/pi*1000*3600^2*24^2;
tau0curl = tau0c*1000*3600^2*24^2*(cos(pi*yp/width)-2/pi);
tau = (taumean+tau0curl)/1000^4/H1; % divide by rho0 and H1 in km - day units
Cr('float','tau',tau);


factor1 = 1000/3600/24*4*1000;
H_0 = -H0*pi/width/factor1;
Hy = H_0*cos(pi*y/width)*pi/width*f0/H1/fr1;

% find h0
if r1 ~= 0
    h0 = H_0*sin(pi*y/width)/r1 + mean(hm(2:end,:),2);

% relaxation conditions
    H_1 = r1*h0; % km /day units
    Hy = diffxy(y,H_1)*f0/H1/fr1;
end
Cr('float','Hy',Hy);

% linear u scaling
uscale = H_0*width/kappa/pi*f0^2/fr1/H1;

if ~exist('q1')
  q=0*(rand(nx,nyp,2)-0.5);
  q(:,1,:)=0;q(:,ny+1,:)=0;
  ubar=zeros(nyp,2);
  ubar(:,1) = -uscale*cos(pi*yp/width)/f0;
  ubar(:,2) = 0;
  
else
    q(:,:,1) = q1';
    q(:,:,2) = q2';
    ubar(:,1) = ubar1;
    ubar(:,2) = ubar2;
end
Cr('float','q',q);



Cr('float','ubar',ubar);



topography = f0/H2*top0*ones(size(yp))*exp(-(x-800).^2/200^2);
%topox = top0*cos(pi*yp/width)*cos(2*x*pi/width)*2*pi/width;
Cr('float','topo',topography');
%Cr('float','topox',topox');

Cr();


% plot barotropic PV contours
%  contourf((f0+beta*yp')*1./(4000-topography(1,:)'/f0*H2*1000))

system('make ch2');
%fid=popen('channelz < channelz.in| tee channelz.out','r')

return

fid=popen('ch2','r')
%fid=fopen('channelphys.out','r');

mycolor;


ts=[];
q1m=[];
q2m=[];
u1m=[];
u2m=[];
qy1m=[];
qy2m=[];
sm=[];
smx=[];

while(1)
  [t,nn]=fread(fid,1,'float');
  if nn ~= 1
    break;
  end
  psi1=fread(fid,[nx,nyp],'float')';
  psi2=fread(fid,[nx,nyp],'float')';
  q1=fread(fid,[nx,nyp],'float')';
  q2=fread(fid,[nx,nyp],'float')';
  u1=fread(fid,[nx,nyp],'float')';
  u2=fread(fid,[nx,nyp],'float')';
  v1=fread(fid,[nx,nyp],'float')';
  v2=fread(fid,[nx,nyp],'float')';
  ubar1=fread(fid,[1,nyp],'float')';
  ubar2=fread(fid,[1,nyp],'float')';
  qbary1=fread(fid,[1,nyp],'float');
  qbary2=fread(fid,[1,nyp],'float');
  phi=fread(fid,[1,nyp],'float')';
  qbar1=cumsum(qbary1)*dx-qbary1(1)*dx/2-qbary1(nyp)/2*dx;
  qbar2=cumsum(qbary2)*dx-qbary2(1)*dx/2-qbary2(nyp)/2*dx;
  qtot1=q1+qbar1'*ones(1,nx);
  qtot2=q2+qbar2'*ones(1,nx);
  psibar1=-cumsum(ubar1)*dx+ubar1(1)*dx/2+ubar1(nyp)/2*dx;
  psibar2=-cumsum(ubar2)*dx+ubar2(1)*dx/2+ubar2(nyp)/2*dx;
  psitot1=psibar1*ones(1,nx)+psi1;
  psitot2=psibar2*ones(1,nx)+psi2;
  %  purge_tmp_files;
  %  title(sprintf('t = %6.2f',t));
  %  contourg(x,yp,qtot1);
  dv=[qtot2;zeros(5,nx);qtot1];
  %  dv=psitot1-psitot2;
  %  dv=[psitot2;zeros(5,nx);psitot1];
  %  dv=[qtot1];
%  df(t,dv,min(min(dv)),max(max(dv)),4);
%  df(t,dv,min(min(dv)),max(max(dv)),1);
  imagesc(dv);axis('xy','equal');title(num2str(t));drawnow
%pause;
  ts=[ts,t];
  q1m=[q1m,[min(min(qtot1));max(max(qtot1))]];
  q2m=[q2m,[min(min(qtot2)),max(max(qtot2))]];
  u1m=[u1m,ubar1];
  u2m=[u2m,ubar2];
  qy1m=[qy1m,qbary1];
  qy2m=[qy2m,qbary2];
  [t,max(ubar1),max(max(q1))]
end
%pclose(fid);
fclose(fid);

q1mm=[min(min(q1m)),max(max(q1m))]
q2mm=[min(min(q2m)),max(max(q2m))]

%plot(ts,max(u1m),ts,max(u2m));
%df(-1);
