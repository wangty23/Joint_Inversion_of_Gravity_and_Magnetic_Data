function qyfw0=fq_fun(A)
[ny,nx]=size(A);
%% 区域数量
qyzx0=[];
nqy=0;
for i=2:ny-1
    for j=2:nx-1
        if(A(i,j)>=A(i-1,j)&&A(i,j)>=A(i,j-1)&&A(i,j)>=A(i+1,j)&&A(i,j)>=A(i,j+1))
            nqy=nqy+1;
            qyzx0=[qyzx0;i,j];
        end
    end
end
%是否有距离很近的�?
qyzx=[];
for i=1:nqy-1
    check=1;
    for j=i+1:nqy
        if(norm(qyzx0(i,:)-qyzx0(j,:))<4)
            check=0;
            break;
        end
    end
    if(check)
        qyzx=[qyzx;qyzx0(i,:)];
    end
end
qyzx=[qyzx;qyzx0(nqy,:)];
%% 各区域所包含的点
nqy=length(qyzx(:,1));
qyfw=zeros(nqy,4);
o=0.3;
for i=1:nqy
    %�?
    for j=qyzx(i,1)-1:-1:2
        if(A(j,qyzx(i,2))<A(j-1,qyzx(i,2))||A(j,qyzx(i,2))<=o*A(qyzx(i,1),qyzx(i,2)))
            qyfw(i,1)=j;
            break;
        end
    end
    %�?
    for j=qyzx(i,1)+1:ny-1
        if(A(j,qyzx(i,2))<A(j+1,qyzx(i,2))||A(j,qyzx(i,2))<=o*A(qyzx(i,1),qyzx(i,2)))
            qyfw(i,2)=j;
            break;
        end
    end
    %�?
    for j=qyzx(i,2)-1:-1:2
        if(A(qyzx(i,1),j)<A(qyzx(i,1),j-1)||A(qyzx(i,1),j)<=o*A(qyzx(i,1),qyzx(i,2)))
            qyfw(i,3)=j;
            break;
        end
    end
    %�?
    for j=qyzx(i,2)+1:nx-1
        if(A(qyzx(i,1),j)<A(qyzx(i,1),j+1)||A(qyzx(i,1),j)<=o*A(qyzx(i,1),qyzx(i,2)))
            qyfw(i,4)=j;
            break;
        end
    end
end
%%
lll=length(qyfw(:,1));
qyfw0=[];
for i=1:lll
    if(abs(qyfw(i,4)-qyfw(i,3))<50&&abs(qyfw(i,2)-qyfw(i,1))<50)
        qyfw0=[qyfw0;qyfw(i,:)];
    end
end
end
















