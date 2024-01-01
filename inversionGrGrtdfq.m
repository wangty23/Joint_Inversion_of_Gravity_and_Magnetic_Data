function [Cg1,Cm1]=inversionGrGrtdfq(Gg,Gt,g,T,Imax,Wdg0,Wdt0,Wzg0,Wzt0,nxs,Dx,Dy,Dz,indqy)
num=0;
k=0;
[aa,bb]=size(Gg);
[~,nqy]=size(indqy);

Wdg1=sparse(1:aa,1:aa,Wdg0); Wdt1=sparse(1:aa,1:aa,Wdt0); 

Wdg1=Wdg1/max(Wdg0); Wdt1=Wdt1/max(Wdt0);
% Wdg1=(Wdg1-min(Wdg0))/(max(Wdg0)-min(Wdg0)); Wdt1=(Wdt1-min(Wdt0))/(max(Wdt0)-min(Wdt0));

% Wdg1=speye(aa); Wdt1=speye(aa); 

Wzg1=sparse(1:bb,1:bb,1./Wzg0); Wzg=sparse(1:bb,1:bb,Wzg0);
Wzt1=sparse(1:bb,1:bb,1./Wzt0); Wzt=sparse(1:bb,1:bb,Wzt0);

Cg1=zeros(bb,1);
Cm1=zeros(bb,1);

g1=Wdg1*g; Vg=Wdg1*Gg*Wzg1; 
T1=Wdt1*T; Vt=Wdt1*Gt*Wzt1; 

rrg=zeros(Imax,1);
rrt=zeros(Imax,1);

ag=1e-1;
at=1e-1;

agg=20;%-1e-3%M对rho
atg=200;%-1e-2%rho对MxMyMz
agg=0;%-1e-3%M对rho
atg=0;
aggx=0;
atgx=0;
aggy=0;
atgy=0;
aggz=0;
atgz=0;

DDx=Dx'*Dx; DDy=Dy'*Dy; DDz=Dz'*Dz;

figure(777)
q111gG=zeros(bb,nqy); q111tG=zeros(bb,nqy);
q111gGx=zeros(bb,nqy); q111gGy=zeros(bb,nqy); q111gGz=zeros(bb,nqy);
q111tGx=zeros(bb,nqy); q111tGy=zeros(bb,nqy); q111tGz=zeros(bb,nqy);
Cgcx=zeros(bb,nqy);    Cgcy=zeros(bb,nqy);    Cgcz=zeros(bb,nqy);
Cmcx=zeros(bb,nqy);    Cmcy=zeros(bb,nqy);    Cmcz=zeros(bb,nqy);
while num<=Imax
    k=k+1;
    disp(100*num/Imax)
    
    if (k>30&&k<80)
        agg=0;%-1e-3%M对rho
        atg=0;
        a1=(1e3);
        a2=5.6*(1e7);
        aggx=a1;
        atgx=a2;
        aggy=a1;
        atgy=a2;
        aggz=a1;
        atgz=a2;
     
    end
    if (k>=80)
         agg=100;%-1e-3%M对rho
        atg=50000;  
    end
    
%        if (mod(k,2))%不能被2整除
%         agg=100;%-1e-3%M对rho
%         atg=50000;
%     else
%         agg=0;%-1e-3%M对rho
%         atg=0;
%         a1=(1e3);
%         a2=4*(1e6);
%         aggx=a1;
%         atgx=a2;
%         aggy=a1;
%         atgy=a2;
%         aggz=a1;
%         atgz=a2;
%        end
    
    Cgc=Cg1; Cmc=Cm1;
    if(k>1)
        Cgc=Cgc/max(abs(Cgc));
        Cmc=Cmc/max(abs(Cmc));
    end
    qgv=0; qtv=0;
    qgvx=0; qgvy=0; qgvz=0; 
    qtvx=0; qtvy=0; qtvz=0;
    for i=1:nqy
        CmCm=(Cmc.*indqy(:,i))'*(Cmc.*indqy(:,i)); 
        CgCg=(Cgc.*indqy(:,i))'*(Cgc.*indqy(:,i)); 
        CmCg=(Cgc.*indqy(:,i))'*(Cmc.*indqy(:,i)); 
        q111gG(:,i)=Wzg*(CmCm*(Cgc.*indqy(:,i))-CmCg*(Cmc.*indqy(:,i)));
        q111tG(:,i)=Wzt*(CgCg*(Cmc.*indqy(:,i))-CmCg*(Cgc.*indqy(:,i)));
        qgv=qgv+q111gG(:,i)'*q111gG(:,i);
        qtv=qtv+q111tG(:,i)'*q111tG(:,i);
        
        Cgci=Cgc.*indqy(:,i); Cmci=Cmc.*indqy(:,i);
        Cgcx(:,i)=Dx*Cgci; Cgcy(:,i)=Dy*Cgci; Cgcz(:,i)=Dz*Cgci;
        Cmcx(:,i)=Dx*Cmci; Cmcy(:,i)=Dy*Cmci; Cmcz(:,i)=Dz*Cmci;
        CmCmx=Cmcx(:,i)'*Cmcx(:,i); CgCgx=Cgcx(:,i)'*Cgcx(:,i); CmCgx=Cmcx(:,i)'*Cgcx(:,i);
        CmCmy=Cmcy(:,i)'*Cmcy(:,i); CgCgy=Cgcy(:,i)'*Cgcy(:,i); CmCgy=Cmcy(:,i)'*Cgcy(:,i);
        CmCmz=Cmcz(:,i)'*Cmcz(:,i); CgCgz=Cgcz(:,i)'*Cgcz(:,i); CmCgz=Cmcz(:,i)'*Cgcz(:,i);
        
        q111gGx(:,i)=Wzg*(CmCmx*(DDx*Cgci)-CmCgx*(DDx*Cmci));
        q111tGx(:,i)=Wzt*(CgCgx*(DDx*Cmci)-CmCgx*(DDx*Cgci));
        q111gGy(:,i)=Wzg*(CmCmy*(DDy*Cgci)-CmCgy*(DDy*Cmci));
        q111tGy(:,i)=Wzt*(CgCgy*(DDy*Cmci)-CmCgy*(DDy*Cgci));
        q111gGz(:,i)=Wzg*(CmCmz*(DDz*Cgci)-CmCgz*(DDz*Cmci));
        q111tGz(:,i)=Wzt*(CgCgz*(DDz*Cmci)-CmCgz*(DDz*Cgci));
        
        qgvx=qgvx+q111gGx(:,i)'*q111gGx(:,i);
        qtvx=qtvx+q111tGx(:,i)'*q111tGx(:,i);
        qgvy=qgvy+q111gGy(:,i)'*q111gGy(:,i);
        qtvy=qtvy+q111tGy(:,i)'*q111tGy(:,i);
        qgvz=qgvz+q111gGz(:,i)'*q111gGz(:,i);
        qtvz=qtvz+q111tGz(:,i)'*q111tGz(:,i);
    end
    
    Cg1=Wzg*Cg1;
    Cm1=Wzt*Cm1;
    
    %     r1g=Vg'*(g1-Vg*Cg1);
    %     r1t=Vt'*(T1-Vt*Cm1);
    %     if(k>1)
    %         r1g=r1g-ag*max(abs(r1g))*Cg1/max(abs(Cg1));
    %         r1t=r1t-at*max(abs(r1t))*Cm1/max(abs(Cm1));
    %     end
    r1g=Vg'*(g1-Vg*Cg1)-ag*Cg1-agg*sum(q111gG,2)-aggx*sum(q111gGx,2)-aggy*sum(q111gGy,2)-aggz*sum(q111gGz,2);
    r1t=Vt'*(T1-Vt*Cm1)-at*Cm1-atg*sum(q111tG,2)-atgx*sum(q111tGx,2)-atgy*sum(q111tGy,2)-atgz*sum(q111tGz,2);
    
    if k==1
        p1g=r1g;
        p1t=r1t;
    else
        u2g=(r1g'*r1g)/(r0g'*r0g);
        p1g=r1g+u2g*p1g;
        u2t=(r1t'*r1t)/(r0t'*r0t);
        p1t=r1t+u2t*p1t;
    end
    r0g=r1g; r0t=r1t;
    %p0g=p1g; p0t=p1t;
    q1g=Vg*p1g; q1t=Vt*p1t;
    q11g=p1g;   q11t=p1t;
    v2g=(r1g'*p1g)/(q1g'*q1g+ag*(q11g'*q11g)+abs(agg)*qgv+abs(aggx)*qgvx+abs(aggy)*qgvy+abs(aggz)*qgvz);
    v2t=(r1t'*p1t)/(q1t'*q1t+at*(q11t'*q11t)+abs(atg)*qtv+abs(atgx)*qtvx+abs(atgy)*qtvy+abs(atgz)*qtvz);
    Cg1=Cg1+v2g*p1g;
    Cm1=Cm1+v2t*p1t;
    Cg1=Wzg1*Cg1;
    Cm1=Wzt1*Cm1;
    Cg1(Cg1>0.5)=0.5; Cg1(Cg1<0)=0;
    Cm1(Cm1>0.5)=0.5; Cm1(Cm1<0)=0;
  
    
    rrg(k,1)=norm(Gg*Cg1-g);
    rrt(k,1)=norm(Gt*Cm1-T);
    subplot(121)
    plot(log10(rrg))
    subplot(122)
    plot(log10(rrt))
    pause(0.001);
    num=num+1;
end
% Cg1=Wzg1*We1_1*Cg1;
% Cm1=Wzt1*We2_1*Cm1;
end