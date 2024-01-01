%% 四面体重力异常正�?
function g=simiantizz(xyz0,xyz)%xyz0:观测点坐�?xyz: 四面体各个角点坐�?
G=6.67*10^-2;
xyzV=zeros(4,6);
xyzV(:,1:3)=[2,3,4;1,3,4;1,2,4;1,2,3];
%% 求外法向
for i=1:4
    p1=xyz(xyzV(i,2),:)-xyz(xyzV(i,1),:);
    p2=xyz(xyzV(i,3),:)-xyz(xyzV(i,1),:);
    p3=xyz(i,:)-xyz(xyzV(i,1),:);
    p4(1)=p1(2)*p2(3)-p1(3)*p2(2);
    p4(2)=p1(3)*p2(1)-p1(1)*p2(3);
    p4(3)=p1(1)*p2(2)-p1(2)*p2(1);
    if(p4*p3'>0)
        xyzV(i,4:6)=-p4/sqrt(p4*p4');
    else
        xyzV(i,4:6)=p4/sqrt(p4*p4');
        ssss=xyzV(i,2);
        xyzV(i,2)=xyzV(i,3);
        xyzV(i,3)=ssss;
    end
end
%% 正演
g=0;
for j=1:4
    %% 夹角
    cosp=xyzV(j,6);
    sinp=sqrt(1-cosp^2);
    if(cosp==-1)
        sint=0;
        cost=1;
    else
        if(cosp==1)
            sint=0;
            cost=-1;
        else
            cost=xyzV(j,4)/sqrt(xyzV(j,4)^2+xyzV(j,5)^2);
            sint=xyzV(j,5)/sqrt(xyzV(j,4)^2+xyzV(j,5)^2);
        end
    end
    if(cosp==0)
        continue;
    end
    %% 1
    for k=1:3
        if(k<3)
            kk=k+1;
        else
            kk=1;
        end
        xc(1)=xyz(xyzV(j,k),1)-xyz0(1);
        xc(2)=xyz(xyzV(j,kk),1)-xyz0(1);
        yc(1)=xyz(xyzV(j,k),2)-xyz0(2);
        yc(2)=xyz(xyzV(j,kk),2)-xyz0(2);
        zc(1)=xyz(xyzV(j,k),3)-xyz0(3);
        zc(2)=xyz(xyzV(j,kk),3)-xyz0(3);
        
        xcc(1)=cost*xc(1)+sint*yc(1);
        xcc(2)=cost*xc(2)+sint*yc(2);
        ycc(1)=cost*yc(1)-sint*xc(1);
        ycc(2)=cost*yc(2)-sint*xc(2);
        zcc(1)=zc(1);
        zcc(2)=zc(2);
        
        xccc(1)=cosp*xcc(1)-sinp*zcc(1);
        xccc(2)=cosp*xcc(2)-sinp*zcc(2);
        yccc(1)=ycc(1);
        yccc(2)=ycc(2);
        zccc(1)=cosp*zcc(1)+sinp*xcc(1);
        zccc(2)=cosp*zcc(2)+sinp*xcc(2);
        
        cosq=(xccc(2)-xccc(1))/sqrt((xccc(2)-xccc(1))^2+(yccc(2)-yccc(1))^2);
        sinq=(yccc(2)-yccc(1))/sqrt((xccc(2)-xccc(1))^2+(yccc(2)-yccc(1))^2);
        
        alpha=-sinp; beta=0; gamma=cosp;
        %alpha=sint*cosp; beta=cost; gamma=-sint*sinp;
        for ii=1:2
            jj=ii;
            uij=(-1)^(ii+1);
            rrrr=sqrt(xccc(ii)^2+yccc(jj)^2+zccc(ii)^2);
            %%
            u=cosq*xccc(ii)+sinq*yccc(jj);
            v=-sinq*xccc(ii)+cosq*yccc(jj);
            w=zccc(ii);
            
%             
            if(u+rrrr<=0)
                a=(alpha*sinq-beta*cosq)*100;
            else
                a=-(alpha*sinq-beta*cosq)*log(u+rrrr);
            end
            
            b=cosq*u*v-(v^2+w^2)*sinq;
            if(w*rrrr*cosq==0)
                if(b<0)
                    c=-(pi/2);
                else
                    c=(pi/2);
                end
            else
                c=atan(b/(w*rrrr*cosq));
            end
            Jji=(a-gamma*c);
            if(isnan(Jji)||isinf(Jji))
                Jji=0;
            end
            %%
            g=g+uij*Jji*cosp;
        end
    end
end
g=-G*g;
end