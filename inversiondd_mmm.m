function [Cg1,Cm1]=inversiondd_mmm(Gg,Gt,g,T,Imax,Wdg0,Wdt0,Wzg0,Wzt0)
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

ag=1e1;
at=1e-2;

figure(777)

while num<=Imax
    k=k+1;
    disp(100*num/Imax)
    
    Cg1=Wzg*Cg1;
    Cm1=Wzt*Cm1;
    
    r1g=Vg'*(g1-Vg*Cg1)-ag*Cg1;
    r1t=Vt'*(T1-Vt*Cm1)-at*Cm1;
    
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
    v2g=(r1g'*p1g)/(q1g'*q1g+ag*(q11g'*q11g));
    v2t=(r1t'*p1t)/(q1t'*q1t+at*(q11t'*q11t));
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