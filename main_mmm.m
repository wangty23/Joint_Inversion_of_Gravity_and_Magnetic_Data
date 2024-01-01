clc
clear
close all
%%
ccc='b_c_d';
load(['Gg' ccc '.mat']);
load(['Gmmm' ccc '.mat']);
load(['Dx' ccc '.mat']);
load(['Dy' ccc '.mat']);
load(['Dz' ccc '.mat']);
xyzg=load(['grav' ccc '.txt'],'-ascii');
xyzm=load(['mag' ccc '.txt'],'-ascii');
Gg=Gg*1000;
%%
g=xyzg(2:end,4);
T=xyzm(2:end,4);
xmin=min(xyzg(2:end,1));
xmax=max(xyzg(2:end,1));
ymin=min(xyzg(2:end,2));
ymax=max(xyzg(2:end,2));
zmin=a;
zmax=0;
nx=b;
ny=c;
nz=d;
N=nx*ny*nz;
Wdm=sum(Gmmm'.^2).^(-0.5);
Wzm=sum(Gmmm.^2).^(0.5);
Wdg=sum(Gg'.^2).^(-0.5);
Wzg=sum(Gg.^2).^(0.5);
Itermax=150;
nxs=2000;
lx=(xmax-xmin)/(nx);
ly=(ymax-ymin)/(ny);
lz=(zmax-zmin)/(nz);
x=(xmin+lx/2):lx:(xmax-lx/2);
y=(ymin+ly/2):ly:(ymax-ly/2);
z=(zmax-lz/2):-lz:(zmin+lz/2);
[xxx,yyy,zzz]=meshgrid(x,y,z);
zzzz=zzz(:);
%%
figure(1111)
nx0=xyzg(1,1); ny0=xyzg(1,2);
gg1=reshape(g,ny0,nx0);
subplot(221)
imagesc(gg1)
colorbar
colormap jet
TT1=reshape(T,ny0,nx0);
subplot(222)
imagesc(TT1)
colorbar
colormap jet
ASS=jxxh(T,nx0,ny0,lx,ly);
subplot(223)
imagesc(ASS)
colorbar
colormap jet
subplot(224)
imagesc(sqrt(abs(ASS).*abs(gg1)))
colorbar
colormap jet
%%
xyz=zeros(N,3);
tt=1;
for i=1:nx
    for j=1:ny
        for k=1:nz
            xyz(tt,:)=[x(i),y(j),z(k)];
            tt=tt+1;
        end
    end
end
Wz0=zeros(1,N);
for tt=1:N
    Wz0(1,tt)=(-xyz(tt,3)-zmin*0.01);
end
Wzg=Wz0.^-1.2;
Wzm0=Wz0.^-1.2;
Wzm=[Wzm0,Wzm0,Wzm0];
%% 
figure(10086)
R=abs(ASS).*abs(gg1);
xx0=reshape(xyzg(2:end,1),ny0,nx0);
yy0=reshape(xyzg(2:end,2),ny0,nx0);
imagesc(xx0(1,:),yy0(:,1),R)
qyfw=fq_fun(R);
[nqy,~]=size(qyfw);
qyfwxy=zeros(size(qyfw));
for i=1:nqy
    qyfwxy(i,:)=[xx0(1,qyfw(i,3)),xx0(1,qyfw(i,4)),yy0(qyfw(i,1),1),yy0(qyfw(i,2),1)];
end
showboxf(qyfwxy)
colorbar
colormap jet
indqy=zeros(N,nqy);
for i=1:nqy
    lim1=xyz(:,1)>=qyfwxy(i,1); lim2=xyz(:,1)<=qyfwxy(i,2);
    lim3=xyz(:,2)>=qyfwxy(i,3); lim4=xyz(:,2)<=qyfwxy(i,4);
    S=(qyfwxy(i,2)-qyfwxy(i,1))*(qyfwxy(i,4)-qyfwxy(i,3));
    indqy(lim1&lim2&lim3&lim4,i)=1;
end
indqy0=ones(N,1)-sum(indqy);
indqy0(indqy0<0)=0;
indqy=[indqy,indqy0];
%%
tic
[Cg,Cm]=inversiondd_mmm(Gg,Gmmm,g,T,Itermax,Wdg,Wdm,Wzg,Wzm);
[Cg,Cm]=inversionGr_mmm(Gg,Gmmm,g,T,Itermax,Wdg,Wdm,Wzg,Wzm);
[Cg,Cm]=inversionGrtd_mmm(Gg,Gmmm,g,T,Itermax,Wdg,Wdm,Wzg,Wzm,Dx,Dy,Dz);
[Cg,Cm]=inversionGr_mmm_fq(Gg,Gmmm,g,T,Itermax,Wdg,Wdm,Wzg,Wzm,indqy);
[Cg,Cm]=inversionGrtd_mmm_fq(Gg,Gmmm,g,T,Itermax,Wdg,Wdm,Wzg,Wzm,Dx,Dy,Dz,indqy);
[Cg,Cm]=inversionGrGrtd_mmm_fq1(Gg,Gmmm,g,T,Itermax,Wdg,Wdm,Wzg,Wzm,Dx,Dy,Dz,indqy);
toc

save('Cg1dd.mat','Cg');
save('Cm1dd.mat','Cm');
%% 
mr=Cg;
r0=Gg*Cg-g;
display(sqrt(r0'*r0/length(g)));
m=Cm;
r0=Gmmm*Cm-T;
display(sqrt(r0'*r0/length(T)));
%%
figure(111)
nx0=xyzg(1,1); ny0=xyzg(1,2);
gg1=reshape(g,ny0,nx0);
subplot(321)
imagesc(gg1)
colorbar
colormap jet
TT1=reshape(T,ny0,nx0);
subplot(322)
imagesc(TT1)
colorbar
colormap jet
gg2=reshape(Gg*Cg,ny0,nx0);
subplot(323)
imagesc(gg2)
colorbar
colormap jet
TT2=reshape(Gmmm*Cm,ny0,nx0);
subplot(324)
imagesc(TT2)
colorbar
colormap jet
subplot(325)
imagesc(gg1-gg2)
colorbar
colormap jet
subplot(326)
imagesc(TT1-TT2)
colorbar
colormap jet
%% mx
Rat=[xmax-xmin,ymax-ymin,zmax-zmin];
load('mod1.mat');
ck=[0,1,1];
mm=reshape(m(1:N,:),nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);

figure(8)
subplot(321)
a=slice(xxx,yyy,zzz,mmm,[],y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(0,0)
%% my
mm=reshape(m(N+1:2*N,:),nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);
subplot(322)
a=slice(xxx,yyy,zzz,mmm,1000,y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(90,0)
%% mz
mm=reshape(m(2*N+1:3*N,:),nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);
subplot(323)
a=slice(xxx,yyy,zzz,mmm,[],y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(0,0)
%% M
subplot(324)
a=slice(xxx,yyy,zzz,mmm,[],y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(0,0)
%% rho
figure(3)
mm=reshape(mr,nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);
a=slice(xxx,yyy,zzz,mmm,[],y(round(ny/2)),[]);
set(a,'EdgeColor','none')
colorbar
colormap jet
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
shading interp
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(0,0)
%%
figure(5)
mm=reshape(sqrt(m(2*N+1:3*N,:).^2+m(N+1:2*N,:).^2+m(1:N,:).^2),nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);
fv1 = isosurface(xxx,yyy,zzz,mmm,0.7*max(sqrt(m(2*N+1:3*N,:).^2+m(N+1:2*N,:).^2+m(1:N,:).^2)));
p1=patch(fv1);
hold on
set(p1,'facecolor','r','edgecolor','none');
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(-20,20)
%%
figure(6)
mm=reshape(mr,nz,ny,nx);
[xxx,yyy,zzz]=meshgrid(x,y,z);
mmm=permute(mm,[2 3 1]);
fv1 = isosurface(xxx,yyy,zzz,mmm,0.7*max(mr));
p1=patch(fv1);
hold on
set(p1,'facecolor','r','edgecolor','none');
axis([min(x),max(x),min(y),max(y),min(z),max(z)])
hold on
grid off
set(gca,'PlotBoxAspectRatio',Rat);
for i=1:length(mod)
    if(mod(i).M==0)
        continue;
    end
    l=length(mod(i).K(:,1));
    for j=1:l
        xk=[mod(i).MD(mod(i).K(j,1),1),mod(i).MD(mod(i).K(j,2),1),mod(i).MD(mod(i).K(j,3),1)];
        yk=[mod(i).MD(mod(i).K(j,1),2),mod(i).MD(mod(i).K(j,2),2),mod(i).MD(mod(i).K(j,3),2)];
        zk=[mod(i).MD(mod(i).K(j,1),3),mod(i).MD(mod(i).K(j,2),3),mod(i).MD(mod(i).K(j,3),3)];
        fill3(xk,yk,zk,ck,'FaceAlpha',0,'LineWidth',1,'EdgeColor',[0 0 0])
    end
    hold on
end
view(-20,20)
%%
figure(7)
plot(mr,sqrt(m(2*N+1:3*N,:).^2+m(N+1:2*N,:).^2+m(1:N,:).^2),'.b')
hold on
plot(0.5,0.5,'r>')
hold on
plot(0.5,0.75,'r^')
hold off