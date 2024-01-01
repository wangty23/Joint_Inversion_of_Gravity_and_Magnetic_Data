clc
clear
close all
%%
xmin=0;
xmax=2000;
ymin=0;
ymax=2000;
nx=b;
ny=c;
lx=(xmax-xmin)/(nx-1);
ly=(ymax-ymin)/(ny-1);
x=xmin:lx:xmax;
y=ymin:ly:ymax;
[xx,yy]=meshgrid(x,y);
zz=0*cos(4*pi*(xx-(xmax+xmin)/2)/(xmax-xmin)).*cos(pi*(yy-(ymax+ymin)/2)/(ymax-ymin));
xyzg=[xx(:),yy(:),zz(:),zz(:)*0];
nxt=b1;
nyt=c1;
lxt=(xmax-xmin)/(nxt-1);
lyt=(ymax-ymin)/(nyt-1);
xt=xmin:lxt:xmax;
yt=ymin:lyt:ymax;
[xxt,yyt]=meshgrid(xt,yt);
zzt=0*cos(4*pi*(xxt-(xmax+xmin)/2)/(xmax-xmin)).*cos(pi*(yyt-(ymax+ymin)/2)/(ymax-ymin));
xyzt=[xxt(:),yyt(:),zzt(:)];
%% Output
outg='grav.txt';
outc='mag.txt';
outtp='topo.txt';
%%
num=2;%the number of models
II0=60;% the inclination of Earth's magnetic field
DD0=45;% the declination of Earth's magnetic field
%%
nnn=1;
mod(nnn).xyz=[550,1000,-260];
mod(nnn).size0=[250,350,200];
t=0;
for i=-1:2:1
    for j=-1:2:1
        for k=-1:2:1
            t=t+1;
             mod(nnn).MD(t,:)=mod(nnn).xyz+[i-0.2*k,j,k].*mod(nnn).size0/2;
        end
    end
end
mod(nnn).rho=0.5;
mod(nnn).M=0.5;
mod(nnn).II=60;
mod(nnn).DD=45;
mod(nnn).s=delaunayTriangulation(mod(nnn).MD);
[mod(nnn).K,~] = convexHull(mod(nnn).s);
%%
nnn=2;
mod(nnn).xyz=[1200,1000,-450];
mod(nnn).size0=[350,400,400];
t=0;
for i=-1:2:1
    for j=-1:2:1
        for k=-1:2:1
            t=t+1;
            mod(nnn).MD(t,:)=mod(nnn).xyz+[i-0.2*i*k,j-0.2*j*k,k].*mod(nnn).size0/2;
        end
    end   
end
mod(nnn).rho=0.5;
mod(nnn).M=0.75;
mod(nnn).II=60;
mod(nnn).DD=45;
mod(nnn).s=delaunayTriangulation(mod(nnn).MD);
[mod(nnn).K,~] = convexHull(mod(nnn).s);
save('mod1.mat','mod');
xyzm=xyzg;
xyzm(:,4)=0;
l=6;
for i=1:num
    Gg=smtK1(xyzg(:,1:3)+1,mod(i).MD,mod(i).s,'Vz');
    Gxx=smtK1(xyzm(:,1:3)+1,mod(i).MD,mod(i).s,'Vxx');
    Gyy=smtK1(xyzm(:,1:3)+1,mod(i).MD,mod(i).s,'Vyy');
    Gzz=smtK1(xyzm(:,1:3)+1,mod(i).MD,mod(i).s,'Vzz');
    Gxz=smtK1(xyzm(:,1:3)+1,mod(i).MD,mod(i).s,'Vxz');
    Gyz=smtK1(xyzm(:,1:3)+1,mod(i).MD,mod(i).s,'Vyz');
    Gxy=smtK1(xyzm(:,1:3)+1,mod(i).MD,mod(i).s,'Vxy');
    Gc=smtKc(Gxx,Gyy,Gzz,Gxz,Gyz,Gxy,mod(i).II,90-mod(i).DD,II0,90-DD0);
    xyzg(:,4)=xyzg(:,4)+Gg*ones(l,1)*mod(i).rho*1000;
    xyzm(:,4)=xyzm(:,4)+Gc*ones(l,1)*mod(i).M;
end
g=reshape(xyzg(:,4),ny,nx);
gm=reshape(xyzm(:,4),ny,nx);

xyzmx=xyzm;xyzmy=xyzm;xyzmz=xyzm;
%%
dep=-1000;
figure(1)
surf(xx,yy,zz,g)
alpha(0.8)
colorbar
colormap jet
hold on
shading interp
ck=[1,0,0];
for i=1:num
    if(mod(i).rho==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',1,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,yk,dep*ones(size(zk)),ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,ymax*ones(size(yk)),zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xmax*ones(size(xk)),yk,zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        for k=1:8
            xxk=mod(i).MD(k,1);yyk=mod(i).MD(k,2);zzk=mod(i).MD(k,3);
            plot3([xxk,xmax],[yyk,yyk],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,ymax],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,yyk],[dep,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
        end
    end
    hold on
end
axis([xmin xmax ymin ymax dep max(max(zz))])
set(gca,'PlotBoxAspectRatio',[2 2 1]);
box on
grid on
%%
figure(777)
surf(xx,yy,zz,gm)
alpha(0.8)
colorbar
colormap jet
hold on
shading interp
ck=[1,0,0];
for i=1:num
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',1,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,yk,dep*ones(size(zk)),ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xk,ymax*ones(size(yk)),zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        fill3(xmax*ones(size(xk)),yk,zk,ck,'FaceAlpha',0,'LineWidth',2,'EdgeColor',[0,0,0])
        for k=1:8
            xxk=mod(i).MD(k,1);yyk=mod(i).MD(k,2);zzk=mod(i).MD(k,3);
            plot3([xxk,xmax],[yyk,yyk],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,ymax],[zzk,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
            plot3([xxk,xxk],[yyk,yyk],[dep,zzk],'--','Color',[0.4,0.4,0.4],'LineWidth',1.5)
        end
    end
    hold on
end
axis([xmin xmax ymin ymax dep max(max(zz))])
set(gca,'PlotBoxAspectRatio',[2 2 1]);
view(-30,20)
box on
grid on
%%
xyzg=[nx,ny,0,0;xyzg];
xyzm=[nx,ny,0,0;xyzm];
xyzt=[nxt,nyt,0;xyzt];
save(outg,'xyzg','-ascii');
save(outc,'xyzm','-ascii');
save(outtp,'xyzt','-ascii');

























