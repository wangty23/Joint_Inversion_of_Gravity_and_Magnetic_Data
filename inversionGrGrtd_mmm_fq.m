function [Cg1,Cm1]=inversionGrGrtd_mmm_fq(Vg,Vt,g,T,Imax,Wdg0,Wdt0,Wzg0,Wzt0,Dx,Dy,Dz,indqy)
num=0;
k=0;
[aa,bb]=size(Vg);
[~,bbt]=size(Vt);
[~,nqy]=size(indqy);

Wdg1=sparse(1:aa,1:aa,Wdg0); Wdt1=sparse(1:aa,1:aa,Wdt0); 

Wdg1=Wdg1/max(Wdg0); Wdt1=Wdt1/max(Wdt0);
% Wdg1=(Wdg1-min(Wdg0))/(max(Wdg0)-min(Wdg0)); Wdt1=(Wdt1-min(Wdt0))/(max(Wdt0)-min(Wdt0));

% Wdg1=speye(aa); Wdt1=speye(aa); 

Wzg1=sparse(1:bb,1:bb,1./Wzg0);   Wzg=sparse(1:bb,1:bb,Wzg0);
Wzt1=sparse(1:bbt,1:bbt,1./Wzt0); Wzt=sparse(1:bbt,1:bbt,Wzt0);

Cg1=zeros(bb,1);
Cm1=zeros(bbt,1);

g1=Wdg1*g; 
Vg=Wdg1*Vg*Wzg1; 
T1=Wdt1*T;
Vt=Wdt1*Vt*Wzt1; 

rrg=zeros(Imax,1);
rrt=zeros(Imax,1);

ag=1e-1;
at=1e-1;


agg=0;  atg=0;  attg=0;

aggx=0; aggy=0; aggz=0;
atgx=0; atgy=0; atgz=0;
attgx=0;attgy=0;attgz=0;

DDx=Dx'*Dx; DDy=Dy'*Dy; DDz=Dz'*Dz;

figure(777)

q111gG=zeros(bb,nqy);  q111tG=zeros(bbt,nqy); q111ttG=zeros(bbt,nqy);

q1gGdx=zeros(bb,nqy);  q1gGdy=zeros(bb,nqy);  q1gGdz=zeros(bb,nqy);
q1tGdx=zeros(bbt,nqy); q1tGdy=zeros(bbt,nqy); q1tGdz=zeros(bbt,nqy);
q1ttGdx=zeros(bbt,nqy);q1ttGdy=zeros(bbt,nqy);q1ttGdz=zeros(bbt,nqy);
    
while num<=Imax
    k=k+1;
    disp(100*num/Imax)
%     if (k>20)
% %         a1=1e3;
% %         a2=1e5;
%         a1=1e3;
%         a2=1e6;
%         aggx=a1;
%         atgx=a2;
%         aggy=a1;
%         atgy=a2;
%         aggz=a1;
%         atgz=a2;
%     end
    if (k>0)
        agg=0.3; atg=0.2; attg=0.2;
        a1=0.01; a2=0.002; a3=0.002;
        a1=0.0; a2=0.00; a3=0.00;
        ax=200;ay=200;az=1;
        aggx=a1; aggy=a1; aggz=a1;
        atgx=ax*a2; atgy=ay*a2; atgz=az*a2;
        attgx=ax*a3;attgy=ay*a3;attgz=az*a3;
    end
    
    Cgc0=Cg1; Cmc0=Cm1;
    if(k>1)
        Cgc0=Cgc0/max(abs(Cgc0));
        Cmc0=Cmc0/max(abs(Cmc0));
    end
    
    qgv=0;   qtv=0;   qttv=0;
    
    qgvdx=0; qgvdy=0; qgvdz=0;
    qtvdx=0; qtvdy=0; qtvdz=0;
    qttvdx=0;qttvdy=0;qttvdz=0;
    
    Cgci=Cgc0;
    Cmcxi=Cmc0(1:bb,1);
    Cmcyi=Cmc0(bb+1:2*bb,1);
    Cmczi=Cmc0(2*bb+1:3*bb,1);
    Cmci=sqrt(Cmcxi.^2+Cmcyi.^2+Cmczi.^2);
    for i=1:nqy
        Cgc=Cgc0.*indqy(:,i); 
        Cmcx=Cmc0(1:bb,1).*indqy(:,i);
        Cmcy=Cmc0(bb+1:2*bb,1).*indqy(:,i);
        Cmcz=Cmc0(2*bb+1:3*bb,1).*indqy(:,i);
        Cmc=sqrt(Cmcx.^2+Cmcy.^2+Cmcz.^2);
        CmCm=Cmc'*Cmc;   CgCg=Cgc'*Cgc;   CmCg=Cmc'*Cgc;
        CmCmx=Cmc'*Cmcx; CmCmy=Cmc'*Cmcy; CmCmz=Cmc'*Cmcz;
        CgCmx=Cgc'*Cmcx; CgCmy=Cgc'*Cmcy; CgCmz=Cgc'*Cmcz;
        q111gG(:,i)=Wzg*(CmCm*Cgc-CmCg*Cmc);
        q111tG(:,i)=Wzt*(CgCg*Cmc0-[CgCmx*Cgc;CgCmy*Cgc;CgCmz*Cgc]);
        q111ttG(:,i)=Wzt*(CmCm*Cmc0-[CmCmx*Cmc;CmCmy*Cmc;CmCmz*Cmc]);
        qgv=qgv+q111gG(:,i)'*q111gG(:,i);
        qtv=qtv+q111tG(:,i)'*q111tG(:,i);
        qttv=qttv+q111ttG(:,i)'*q111ttG(:,i);
        
        Cgcdx=(Dx*Cgci).*indqy(:,i);   Cgcdy=(Dy*Cgci).*indqy(:,i);   Cgcdz=(Dz*Cgci).*indqy(:,i);%密度梯度
        Cmcdx=(Dx*Cmci).*indqy(:,i);   Cmcdy=(Dy*Cmci).*indqy(:,i);   Cmcdz=(Dz*Cmci).*indqy(:,i);%M梯度
        Cmcxdx=(Dx*Cmcxi).*indqy(:,i); Cmcxdy=(Dy*Cmcxi).*indqy(:,i); Cmcxdz=(Dz*Cmcxi).*indqy(:,i);%Mx梯度
        Cmcydx=(Dx*Cmcyi).*indqy(:,i); Cmcydy=(Dy*Cmcyi).*indqy(:,i); Cmcydz=(Dz*Cmcyi).*indqy(:,i);%My梯度
        Cmczdx=(Dx*Cmczi).*indqy(:,i); Cmczdy=(Dy*Cmczi).*indqy(:,i); Cmczdz=(Dz*Cmczi).*indqy(:,i);%Mz梯度
        
        CgCgdx=Cgcdx'*Cgcdx;   CgCgdy=Cgcdy'*Cgcdy;   CgCgdz=Cgcdz'*Cgcdz;
        CmCmdx=Cmcdx'*Cmcdx;   CmCmdy=Cmcdy'*Cmcdy;   CmCmdz=Cmcdz'*Cmcdz;
        CmCgdx=Cmcdx'*Cgcdx;   CmCgdy=Cmcdy'*Cgcdy;   CmCgdz=Cmcdz'*Cgcdz;
        CmCmxdx=Cmcdx'*Cmcxdx; CmCmxdy=Cmcdy'*Cmcxdy; CmCmxdz=Cmcdz'*Cmcxdz;
        CmCmydx=Cmcdx'*Cmcydx; CmCmydy=Cmcdy'*Cmcydy; CmCmydz=Cmcdz'*Cmcydz;
        CmCmzdx=Cmcdx'*Cmczdx; CmCmzdy=Cmcdy'*Cmczdy; CmCmzdz=Cmcdz'*Cmczdz;
        CgCmxdx=Cgcdx'*Cmcxdx; CgCmxdy=Cgcdy'*Cmcxdy; CgCmxdz=Cgcdz'*Cmcxdz;
        CgCmydx=Cgcdx'*Cmcydx; CgCmydy=Cgcdy'*Cmcydy; CgCmydz=Cgcdz'*Cmcydz;
        CgCmzdx=Cgcdx'*Cmczdx; CgCmzdy=Cgcdy'*Cmczdy; CgCmzdz=Cgcdz'*Cmczdz;
        
        %重力
        DCgcdx=(DDx*Cgci).*indqy(:,i); DCgcdy=(DDy*Cgci).*indqy(:,i); DCgcdz=(DDz*Cgci).*indqy(:,i);
        DCmcdx=(DDx*Cmci).*indqy(:,i); DCmcdy=(DDy*Cmci).*indqy(:,i); DCmcdz=(DDz*Cmci).*indqy(:,i);
        q1gGdx(:,i)=Wzg*(CmCmdx*DCgcdx-CmCgdx*DCmcdx);
        q1gGdy(:,i)=Wzg*(CmCmdy*DCgcdy-CmCgdy*DCmcdy);
        q1gGdz(:,i)=Wzg*(CmCmdz*DCgcdz-CmCgdz*DCmcdz);
        %rho->MxMyMz
        qqqqx=[(DDx*Cmcxi).*indqy(:,i);(DDx*Cmcyi).*indqy(:,i);(DDx*Cmczi).*indqy(:,i)];
        qqqqy=[(DDy*Cmcxi).*indqy(:,i);(DDy*Cmcyi).*indqy(:,i);(DDy*Cmczi).*indqy(:,i)];
        qqqqz=[(DDz*Cmcxi).*indqy(:,i);(DDz*Cmcyi).*indqy(:,i);(DDz*Cmczi).*indqy(:,i)];
        q1tGdx(:,i)=Wzt*(CgCgdx*qqqqx...
            -[CgCmxdx*DCgcdx;CgCmydx*DCgcdx;CgCmzdx*DCgcdx]);
        q1tGdy(:,i)=Wzt*(CgCgdy*qqqqy...
            -[CgCmxdy*DCgcdy;CgCmydy*DCgcdy;CgCmzdy*DCgcdy]);
        q1tGdz(:,i)=Wzt*(CgCgdz*qqqqz...
            -[CgCmxdz*DCgcdz;CgCmydz*DCgcdz;CgCmzdz*DCgcdz]);
        %M->MxMyMz
        q1ttGdx(:,i)=Wzt*(CmCmdx*qqqqx...
            -[CmCmxdx*DCmcdx;CmCmydx*DCmcdx;CmCmzdx*DCmcdx]);
        q1ttGdy(:,i)=Wzt*(CmCmdy*qqqqy...
            -[CmCmxdy*DCmcdy;CmCmydy*DCmcdy;CmCmzdy*DCmcdy]);
        q1ttGdz(:,i)=Wzt*(CmCmdz*qqqqz...
            -[CmCmxdz*DCgcdz;CmCmydz*DCmcdz;CmCmzdz*DCgcdz]);
        
        %
        qgvdx=qgvdx+q1gGdx(:,i)'*q1gGdx(:,i);
        qgvdy=qgvdy+q1gGdy(:,i)'*q1gGdy(:,i);
        qgvdz=qgvdz+q1gGdz(:,i)'*q1gGdz(:,i);
        qtvdx=qtvdx+q1tGdx(:,i)'*q1tGdx(:,i);
        qtvdy=qtvdy+q1tGdy(:,i)'*q1tGdy(:,i);
        qtvdz=qtvdz+q1tGdz(:,i)'*q1tGdz(:,i);
        qttvdx=qttvdx+q1ttGdx(:,i)'*q1ttGdx(:,i);
        qttvdy=qttvdy+q1ttGdy(:,i)'*q1ttGdy(:,i);
        qttvdz=qttvdz+q1ttGdz(:,i)'*q1ttGdz(:,i);
    end
    
    Cg1=Wzg*Cg1;
    Cm1=Wzt*Cm1;
    
    if(k>1)
        agg1=agg*norm(r1g0)/norm(sum(q111gG,2));
        atg1=atg*norm(r1t0)/norm(sum(q111tG,2));
        attg1=attg*norm(r1t0)/norm(sum(q111ttG,2));
        
        aggx1=aggx*norm(r1g0)/norm(sum(q1gGdx,2));  
        aggy1=aggy*norm(r1g0)/norm(sum(q1gGdy,2));  
        aggz1=aggz*norm(r1g0)/norm(sum(q1gGdz,2));
        atgx1=atgx*norm(r1t0)/(norm(sum(q1tGdx,2))+norm(sum(q1tGdy,2))+norm(sum(q1tGdz,2)));
        atgy1=atgy*norm(r1t0)/(norm(sum(q1tGdx,2))+norm(sum(q1tGdy,2))+norm(sum(q1tGdz,2)));
        atgz1=atgz*norm(r1t0)/(norm(sum(q1tGdx,2))+norm(sum(q1tGdy,2))+norm(sum(q1tGdz,2)));
        attgx1=attgx*norm(r1t0)/(norm(sum(q1ttGdx,2))+norm(sum(q1ttGdy,2))+norm(sum(q1ttGdz,2)));
        attgy1=attgy*norm(r1t0)/(norm(sum(q1ttGdx,2))+norm(sum(q1ttGdy,2))+norm(sum(q1ttGdz,2)));
        attgz1=attgz*norm(r1t0)/(norm(sum(q1ttGdx,2))+norm(sum(q1ttGdy,2))+norm(sum(q1ttGdz,2)));
    else
        agg1=agg;    atg1=atg;    attg1=attg;

        aggx1=aggx;  aggy1=aggy;  aggz1=aggz;
        atgx1=atgx;  atgy1=atgy;  atgz1=atgz;
        attgx1=attgx;attgy1=attgy;attgz1=attgz;
    end
    r1g0=Vg'*(g1-Vg*Cg1)-ag*Cg1;
    r1t0=Vt'*(T1-Vt*Cm1)-at*Cm1;
    
    r1g=r1g0-agg1*sum(q111gG,2)-aggx1*sum(q1gGdx,2)-aggy1*sum(q1gGdy,2)-aggz1*sum(q1gGdz,2);
    r1t=r1t0-atg1*sum(q111tG,2)-attg1*sum(q111ttG,2)...
        -atgx1*sum(q1tGdx,2)-atgy1*sum(q1tGdy,2)-atgz1*sum(q1tGdz,2)...
        -attgx1*sum(q1ttGdx,2)-attgy1*sum(q1ttGdy,2)-attgz1*sum(q1ttGdz,2);
    
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
    qqg=abs(aggx1)*qgvdx+abs(aggy1)*qgvdy+abs(aggz1)*qgvdz;
    qqt=abs(atgx1)*qtvdx+abs(atgy1)*qtvdy+abs(atgz1)*qtvdz...
        +abs(attgx1)*qttvdx+abs(attgy1)*qttvdy+abs(attgz1)*qttvdz;
    v2g=(r1g'*p1g)/(q1g'*q1g+ag*(q11g'*q11g)+abs(agg1)*qgv+qqg);
    v2t=(r1t'*p1t)/(q1t'*q1t+at*(q11t'*q11t)+abs(atg1)*qtv+abs(attg1)*qttv+qqt);
    Cg1=Cg1+v2g*p1g;
    Cm1=Cm1+v2t*p1t;
    Cg1=Wzg1*Cg1;
    Cm1=Wzt1*Cm1;
    Cg1(Cg1>0.5)=0.5; Cg1(Cg1<0)=0;
    Cm1(Cm1>0.5)=0.5; Cm1(Cm1<-0.5)=-0.5;
    
    rrg(k,1)=norm(r0g);
    rrt(k,1)=norm(r0t);
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