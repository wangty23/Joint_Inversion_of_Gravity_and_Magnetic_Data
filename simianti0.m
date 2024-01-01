%% 四面体重力异常正�?
function g=simianti0(xyz0,xyz,type)%xyz0:观测点坐�?xyz: 四面体各个角点坐�?
switch(type)
    case 'Vz'
        g=simianti(xyz0,xyz);
    case 'Vx'
        AA=[0,0,-1;0,1,0;1,0,0];%x
        xyz0=xyz0*AA;
        xyz=xyz*AA;
        g=simianti(xyz0,xyz);
    case 'Vy'
        AA=[1,0,0;0,0,-1;0,1,0];%y
        xyz0=xyz0*AA;
        xyz=xyz*AA;
        g=simianti(xyz0,xyz);
    case 'Vzz'
         g=simiantizz(xyz0,xyz);
         % g=Gzz(xyz0,xyz);
    case 'Vxx'
        AA=[0,0,-1;0,1,0;1,0,0];%x
        xyz0=xyz0*AA;
        xyz=xyz*AA;
        g=simiantizz(xyz0,xyz);
    case 'Vyy'
        AA=[1,0,0;0,0,-1;0,1,0];%y
        xyz0=xyz0*AA;
        xyz=xyz*AA;
        g=simiantizz(xyz0,xyz);
    case 'Vyz'
        g=simiantiyz(xyz0,xyz);
    case 'Vxy'
        AA=[0,0,-1;0,1,0;1,0,0];%x
        xyz0=xyz0*AA;
        xyz=xyz*AA;
        g=simiantiyz(xyz0,xyz);
    case 'Vxz'
        AA=[0,1,0;1,0,0;0,0,1];%y
        xyz0=xyz0*AA;
        xyz=xyz*AA;
        g=simiantiyz(xyz0,xyz);
end
end