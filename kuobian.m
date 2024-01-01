%% 余弦扩边函数
function [ddkb,Nkb,NxL,NyL,NxR,NyR]=kuobian(dd,nx,ny)
Nkb=2^(floor(log2(max([nx,ny])))+1);%计算扩边点数
NxL=floor((Nkb-nx)/2);
NyL=floor((Nkb-ny)/2);
NxR=Nkb-nx-NxL;
NyR=Nkb-ny-NyL;
ddkb=zeros(Nkb);
ddkb(NyL+1:NyL+ny,NxL+1:NxL+nx)=dd;
%左
for i=NyL+1:NyL+ny
    for j=1:NxL
        ddkb(i,j)=ddkb(i,NxL+1)*(1+cos(pi*(NxL+1-j)/NxL))/2;
    end
end
%右
for i=NyL+1:NyL+ny
    for j=1:NxR
        ddkb(i,NxL+nx+j)=ddkb(i,NxL+nx)*(1+cos(pi*(j)/NxR))/2;
    end
end
%下
for i=NxL+1:NxL+nx
    for j=1:NyL
        ddkb(j,i)=ddkb(NyL+1,i)*(1+cos(pi*(NyL+1-j)/NyL))/2;
    end
end
%上
for i=NxL+1:NxL+nx
    for j=1:NyR
        ddkb(NyL+ny+j,i)=ddkb(NyL+ny,i)*(1+cos(pi*(j)/NyR))/2;
    end
end
end