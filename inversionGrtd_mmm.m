function [Cg1,Cm1]=inversionGrtd_mmm(Gg,Gt,g,T,Imax,Wdg0,Wdt0,Wzg0,Wzt0,Dx,Dy,Dz)
num=0;
k=0;
[aa,bb]=size(Gg);
[~,bbt]=size(Gt);

Wdg1=sparse(1:aa,1:aa,Wdg0); Wdt1=sparse(1:aa,1:aa,Wdt0); 

Wdg1=Wdg1/max(Wdg0); Wdt1=Wdt1/max(Wdt0);


Wzg1=sparse(1:bb,1:bb,1./Wzg0);   Wzg=sparse(1:bb,1:bb,Wzg0);
Wzt1=sparse(1:bbt,1:bbt,1./Wzt0); Wzt=sparse(1:bbt,1:bbt,Wzt0);

Cg1=zeros(bb,1);
Cm1=zeros(bbt,1);

g1=Wdg1*g; Vg=Wdg1*Gg*Wzg1; 
T1=Wdt1*T; Vt=Wdt1*Gt*Wzt1; 

rrg=zeros(Imax,1);
rrt=zeros(Imax,1);

ag=1e-1;
at=1e-1;


aggx=0; aggy=0; aggz=0;
atgx=0; atgy=0; atgz=0;
attgx=0;attgy=0;attgz=0;

figure(777)

while num<=Imax
    k=k+1;
    disp(100*num/Imax)

    if (k>0)
        a1=0.8; a2=0.45; a3=0.45;
        aggx=a1; aggy=a1; aggz=a1;
        atgx=a2; atgy=a2; atgz=a2;
        attgx=a3;attgy=a3;attgz=a3;
    end
    
    Cgc=Cg1; Cmc0=Cm1;
    if(k>1)
        Cgc=Cgc/max(abs(Cgc));
        Cmc0=Cmc0/max(abs(Cmc0));
    end
    Cmcx=Cmc0(1:bb,1);        Cmcy=Cmc0(bb+1:2*bb,1); 
    Cmcz=Cmc0(2*bb+1:3*bb,1); Cmc=sqrt(Cmcx.^2+Cmcy.^2+Cmcz.^2);
    
    Cgcdx=Dx*Cgc;   Cgcdy=Dy*Cgc;   Cgcdz=Dz*Cgc;
    Cmcdx=Dx*Cmc;   Cmcdy=Dy*Cmc;   Cmcdz=Dz*Cmc;
    Cmcxdx=Dx*Cmcx; Cmcxdy=Dy*Cmcx; Cmcxdz=Dz*Cmcx;
    Cmcydx=Dx*Cmcy; Cmcydy=Dy*Cmcy; Cmcydz=Dz*Cmcy;
    Cmczdx=Dx*Cmcz; Cmczdy=Dy*Cmcz; Cmczdz=Dz*Cmcz;
    
    CgCgdx=Cgcdx'*Cgcdx;   CgCgdy=Cgcdy'*Cgcdy;   CgCgdz=Cgcdz'*Cgcdz;
    CmCmdx=Cmcdx'*Cmcdx;   CmCmdy=Cmcdy'*Cmcdy;   CmCmdz=Cmcdz'*Cmcdz;
    CmCgdx=Cmcdx'*Cgcdx;   CmCgdy=Cmcdy'*Cgcdy;   CmCgdz=Cmcdz'*Cgcdz;
    CmCmxdx=Cmcdx'*Cmcxdx; CmCmxdy=Cmcdy'*Cmcxdy; CmCmxdz=Cmcdz'*Cmcxdz;
    CmCmydx=Cmcdx'*Cmcydx; CmCmydy=Cmcdy'*Cmcydy; CmCmydz=Cmcdz'*Cmcydz;
    CmCmzdx=Cmcdx'*Cmczdx; CmCmzdy=Cmcdy'*Cmczdy; CmCmzdz=Cmcdz'*Cmczdz;
    CgCmxdx=Cgcdx'*Cmcxdx; CgCmxdy=Cgcdy'*Cmcxdy; CgCmxdz=Cgcdz'*Cmcxdz;
    CgCmydx=Cgcdx'*Cmcydx; CgCmydy=Cgcdy'*Cmcydy; CgCmydz=Cgcdz'*Cmcydz;
    CgCmzdx=Cgcdx'*Cmczdx; CgCmzdy=Cgcdy'*Cmczdy; CgCmzdz=Cgcdz'*Cmczdz;
    
    
    q1gGdx=Wzg*(CmCmdx*(Dx'*Cgcdx)-CmCgdx*(Dx'*Cmcdx)); 
    q1gGdy=Wzg*(CmCmdy*(Dy'*Cgcdy)-CmCgdy*(Dy'*Cmcdy)); 
    q1gGdz=Wzg*(CmCmdz*(Dz'*Cgcdz)-CmCgdz*(Dz'*Cmcdz)); 
    %rho->MxMyMz
    DCgcdx=Dx'*Cgcdx; DCgcdy=Dy'*Cgcdy; DCgcdz=Dz'*Cgcdz;
    q1tGdx=Wzt*(CgCgdx*[Dx'*Cmcxdx;Dx'*Cmcydx;Dx'*Cmczdx]...
        -[CgCmxdx*DCgcdx;CgCmydx*DCgcdx;CgCmzdx*DCgcdx]);
    q1tGdy=Wzt*(CgCgdy*[Dy'*Cmcxdy;Dy'*Cmcydy;Dy'*Cmczdy]...
        -[CgCmxdy*DCgcdy;CgCmydy*DCgcdy;CgCmzdy*DCgcdy]);
    q1tGdz=Wzt*(CgCgdz*[Dz'*Cmcxdz;Dz'*Cmcydz;Dz'*Cmczdz]...
        -[CgCmxdz*DCgcdz;CgCmydz*DCgcdz;CgCmzdz*DCgcdz]);
    %M->MxMyMz
    DCmcdx=Dx'*Cmcdx; DCmcdy=Dy'*Cmcdy; DCmcdz=Dz'*Cmcdz;
    q1ttGdx=Wzt*(CmCmdx*[Dx'*Cmcxdx;Dx'*Cmcydx;Dx'*Cmczdx]...
        -[CmCmxdx*DCmcdx;CmCmydx*DCmcdx;CmCmzdx*DCmcdx]);
    q1ttGdy=Wzt*(CmCmdy*[Dy'*Cmcxdy;Dy'*Cmcydy;Dy'*Cmczdy]...
        -[CmCmxdy*DCmcdy;CmCmydy*DCmcdy;CmCmzdy*DCmcdy]);
    q1ttGdz=Wzt*(CmCmdz*[Dz'*Cmcxdz;Dz'*Cmcydz;Dz'*Cmczdz]...
        -[CmCmxdz*DCgcdz;CmCmydz*DCmcdz;CmCmzdz*DCgcdz]);
    
    Cg1=Wzg*Cg1;
    Cm1=Wzt*Cm1;
    
    if(k>1)
        aggx1=aggx*norm(r1g0)/norm(q1gGdx);  
        aggy1=aggy*norm(r1g0)/norm(q1gGdy);  
        aggz1=aggz*norm(r1g0)/norm(q1gGdz);
        atgx1=atgx*norm(r1t0)/(norm(q1tGdx)+norm(q1tGdy)+norm(q1tGdz));
        atgy1=atgy*norm(r1t0)/(norm(q1tGdx)+norm(q1tGdy)+norm(q1tGdz));
        atgz1=atgz*norm(r1t0)/(norm(q1tGdx)+norm(q1tGdy)+norm(q1tGdz));
        attgx1=attgx*norm(r1t0)/(norm(q1ttGdx)+norm(q1ttGdy)+norm(q1ttGdz));
        attgy1=attgy*norm(r1t0)/(norm(q1ttGdx)+norm(q1ttGdy)+norm(q1ttGdz));
        attgz1=attgz*norm(r1t0)/(norm(q1ttGdx)+norm(q1ttGdy)+norm(q1ttGdz));
    else
        aggx1=aggx;  aggy1=aggy;  aggz1=aggz;
        atgx1=atgx;  atgy1=atgy;  atgz1=atgz;
        attgx1=attgx;attgy1=attgy;attgz1=attgz;
    end
    r1g0=Vg'*(g1-Vg*Cg1)-ag*Cg1;
    r1t0=Vt'*(T1-Vt*Cm1)-at*Cm1;
    
    r1g=r1g0-aggx1*q1gGdx-aggy1*q1gGdy-aggz1*q1gGdz;
    r1t=r1t0-atgx1*q1tGdx-atgy1*q1tGdy-atgz1*q1tGdz...
        -attgx1*q1ttGdx-attgy1*q1ttGdy-attgz1*q1ttGdz;
    
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
    q1g=Vg*p1g; q1t=Vt*p1t;
    q11g=p1g;   q11t=p1t;
    qqg=abs(aggx1)*(q1gGdx'*q1gGdx)+abs(aggy1)*(q1gGdy'*q1gGdy)+abs(aggz1)*(q1gGdz'*q1gGdz);
    qqt=abs(atgx1)*(q1tGdx'*q1tGdx)+abs(atgy1)*(q1tGdy'*q1tGdy)+abs(atgz1)*(q1tGdz'*q1tGdz)...
        +abs(attgx1)*(q1ttGdx'*q1ttGdx)+abs(attgy1)*(q1ttGdy'*q1ttGdy)+abs(attgz1)*(q1ttGdz'*q1ttGdz);
    v2g=(r1g'*p1g)/(q1g'*q1g+ag*(q11g'*q11g)+qqg);
    v2t=(r1t'*p1t)/(q1t'*q1t+at*(q11t'*q11t)+qqt);
    Cg1=Cg1+v2g*p1g;
    Cm1=Cm1+v2t*p1t;
    Cg1=Wzg1*Cg1;
    Cm1=Wzt1*Cm1;
    Cg1(Cg1>0.5)=0.5; Cg1(Cg1<0)=0;
    Cm1(Cm1>0.5)=0.5; Cm1(Cm1<-0.5)=-0.5;
    
    rrg(k,1)=norm(Gg*Cg1-g);
    rrt(k,1)=norm(Gt*Cm1-T);
    subplot(121)
    plot(log10(rrg))
    subplot(122)
    plot(log10(rrt))
    pause(0.001);
    num=num+1;
end

end