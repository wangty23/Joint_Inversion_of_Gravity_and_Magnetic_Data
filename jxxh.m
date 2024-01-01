function ASS=jxxh(T,nx,ny,lx,ly)
dd=reshape(T,ny,nx);%数据数组恢复成nx×ny的网格数据
[ddkb,Nkb,NxL,NyL,NxR,NyR]=kuobian(dd,nx,ny);
%% 波数
u=fftshift((1:Nkb)-Nkb/2-1);%x方向波数
u=-2*pi.*(u)./((Nkb-1)*lx);
v=fftshift((1:Nkb)-Nkb/2-1);%y方向波数
v=-2*pi.*(v)./((Nkb-1)*ly);
[uu,vv]=meshgrid(u,v);
w=sqrt(uu.^2+vv.^2);
%% FFT
Fddkb=fft2(ddkb);
%% 求导
[ddkbx,ddkby]=gradient(ddkb);
%x
ddx=ddkbx(NyL+1:NyL+ny,NxL+1:NxL+nx)/lx;
%y
ddy=ddkby(NyL+1:NyL+ny,NxL+1:NxL+nx)/ly;
%z
Fddkbz=w.*Fddkb;
ddkbz=real(ifft2(Fddkbz));
ddz=ddkbz(NyL+1:NyL+ny,NxL+1:NxL+nx);
%%
ASS=sqrt(ddx.^2+ddy.^2+ddz.^2);
end