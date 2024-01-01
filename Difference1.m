function [Dx1,Dy1,Dz1]=Difference1(nx,ny,nz)
N=nx*ny*nz;
Dx=zeros(2*N,3);
Dy=zeros(2*N,3);
Dz=zeros(2*N,3);
for i=0:nx-1
    for j=0:ny-1
        for k=0:nz-1
            %% z
            t=k+1+j*nz+i*ny*nz;
            if(k==0)
                Dz(t+0,:)=[t,t+0,-1];
                Dz(t+N,:)=[t,t+1,+1];
            else
                if(k==nz-1)
                    Dz(t+0,:)=[t,t+0,+1];
                    Dz(t+N,:)=[t,t-1,-1];
                else
                    Dz(t+0,:)=[t,t+1,+1/2];
                    Dz(t+N,:)=[t,t-1,-1/2];
                end
            end
            %% y
            if(j==0)
                Dy(t+0,:)=[t,t+00,-1];
                Dy(t+N,:)=[t,t+nz,+1];
            else
                if(j==ny-1)
                    Dy(t+0,:)=[t,t+00,+1];
                    Dy(t+N,:)=[t,t-nz,-1];
                else
                    Dy(t+0,:)=[t,t+nz,+1/2];
                    Dy(t+N,:)=[t,t-nz,-1/2];
                end
            end
            %% x
            if(i==0)
                Dx(t+0,:)=[t,t+00000,-1];
                Dx(t+N,:)=[t,t+ny*nz,+1];
            else
                if(i==nx-1)
                    Dx(t+0,:)=[t,t+00000,+1];
                    Dx(t+N,:)=[t,t-ny*nz,-1];
                else
                    Dx(t+0,:)=[t,t+ny*nz,+1/2];
                    Dx(t+N,:)=[t,t-ny*nz,-1/2];
                end
            end
        end
    end
end
Dx1=sparse(Dx(:,1),Dx(:,2),Dx(:,3));
Dy1=sparse(Dy(:,1),Dy(:,2),Dy(:,3));
Dz1=sparse(Dz(:,1),Dz(:,2),Dz(:,3));
end