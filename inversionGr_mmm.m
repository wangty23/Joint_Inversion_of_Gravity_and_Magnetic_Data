function [Cg1,Cm1]=inversionGr_mmm(Gg,Gt,g,T,Imax,Wdg0,Wdt0,Wzg0,Wzt0)
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

agg=0;
atg=0;
attg=0;

figure(777)

while num<=Imax
    k=k+1;
    disp(100*num/Imax)
    if (k>0)
        agg=0.5;
        atg=0.25;
        attg=0.25;
    end
    
    Cgc=Cg1; Cmc0=Cm1;
    if(k>1)
        Cgc=Cgc/max(abs(Cgc));
        Cmc0=Cmc0/max(abs(Cmc0));
    end
    Cmcx=Cmc0(1:bb,1);        Cmcy=Cmc0(bb+1:2*bb,1); 
    Cmcz=Cmc0(2*bb+1:3*bb,1); Cmc=sqrt(Cmcx.^2+Cmcy.^2+Cmcz.^2);
    CmCm=Cmc'*Cmc;   CgCg=Cgc'*Cgc;   CmCg=Cmc'*Cgc;
    CmCmx=Cmc'*Cmcx; CmCmy=Cmc'*Cmcy; CmCmz=Cmc'*Cmcz;
    CgCmx=Cgc'*Cmcx; CgCmy=Cgc'*Cmcy; CgCmz=Cgc'*Cmcz;
    q111gG=Wzg*(CmCm*Cgc-CmCg*Cmc);
    q111tG=Wzt*(CgCg*Cmc0-[CgCmx*Cgc;CgCmy*Cgc;CgCmz*Cgc]);
    q111ttG=Wzt*(CmCm*Cmc0-[CmCmx*Cmc;CmCmy*Cmc;CmCmz*Cmc]);
    
    Cg1=Wzg*Cg1;
    Cm1=Wzt*Cm1;
    
    
    if(k>1)
        agg1=agg*norm(r1g0)/norm(q111gG);
        atg1=atg*norm(r1t0)/norm(q111tG);
        attg1=attg*norm(r1t0)/norm(q111ttG);
    else
        agg1=agg;
        atg1=atg;
        attg1=attg;
    end
    r1g0=Vg'*(g1-Vg*Cg1)-ag*Cg1;
    r1t0=Vt'*(T1-Vt*Cm1)-at*Cm1;
    r1g=r1g0-agg1*q111gG;
    r1t=r1t0-atg1*q111tG-attg1*q111ttG;
    
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
    v2g=(r1g'*p1g)/(q1g'*q1g+ag*(q11g'*q11g)+abs(agg1)*(q111gG'*q111gG));
    v2t=(r1t'*p1t)/(q1t'*q1t+at*(q11t'*q11t)+abs(atg1)*(q111tG'*q111tG)+abs(attg1)*(q111ttG'*q111ttG));
    Cg1=Cg1+v2g*p1g;
    Cm1=Cm1+v2t*p1t;
    Cg1=Wzg1*Cg1;
    Cm1=Wzt1*Cm1;
    
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