function G=smtK1(xyz,MD,MDp,type)
nobs=length(xyz(:,1));
N=length(MDp(:,1));
G=zeros(nobs,N);
xyz_cor=zeros(4,3,N);
for j=1:N
    xyz_cor(:,:,j)=MD(MDp(j,:),:);
end
for i=1:nobs
    disp(i);
    xyz_obs=xyz(i,:);
    %% 
    for j=1:N
        G(i,j)=simianti0(xyz_obs,xyz_cor(:,:,j),type);
    end
end
end
