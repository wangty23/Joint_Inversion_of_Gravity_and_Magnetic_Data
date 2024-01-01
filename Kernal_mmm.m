clc
clear
close all
ccc='b_c_d';
xyzg=load(['grav' ccc '.txt'],'-ascii');
outg=['Gg' ccc '.mat'];
% outm=['Gm' ccc '.mat'];
outmmm=['Gmmm' ccc '.mat'];
outDx=['Dx' ccc '.mat'];
outDy=['Dy' ccc '.mat'];
outDz=['Dz' ccc '.mat'];
nobs=length(xyzg(:,1))-1;
xmin=min(xyzg(2:end,1));
xmax=max(xyzg(2:end,1));
ymin=min(xyzg(2:end,2));
ymax=max(xyzg(2:end,2));
zmin=a ;
zmax=0;
nx=b;
ny=c;
nz=d;
N=nx*ny*nz;
II0=e;
DD0=f;
II=g;
DD=h;
lx=(xmax-xmin)/(nx);
ly=(ymax-ymin)/(ny);
lz=(zmax-zmin)/(nz);
size0=[lx,ly,lz];
x=(xmin+lx/2):lx:(xmax-lx/2);
y=(ymin+ly/2):ly:(ymax-ly/2);
z=(zmax-lz/2):-lz:(zmin+lz/2);
Gg=zeros(nobs,N);
Gm=zeros(nobs,N); Gmx=zeros(nobs,N);
Gmy=zeros(nobs,N);Gmz=zeros(nobs,N);
tt=1;
xyz=zeros(N,3);
for i=1:nx
    for j=1:ny
        for k=1:nz
            xyz(tt,:)=[x(i),y(j),z(k)];
            tt=tt+1;
        end
    end
end
for t=1:nobs
    xyz0=xyzg(t+1,1:3);
    parfor tt=1:N
        Gg(t,tt)=grav_fun('Vz',xyz0,xyz(tt,:),size0);
        [Gm(t,tt),Gmx(t,tt),Gmy(t,tt),Gmz(t,tt)]=mag_fun_mmm('Vt',xyz0,xyz(tt,:),size0,II,90-DD,II0,90-DD0);%磁mmm
    end
    disp(t);
end
Gmmm=[Gmx,Gmy,Gmz];
save(outg,'Gg');
save(outmmm,'Gmmm');
[Dx,Dy,Dz]=Difference1(nx,ny,nz);
save(outDx,'Dx');
save(outDy,'Dy');
save(outDz,'Dz');












