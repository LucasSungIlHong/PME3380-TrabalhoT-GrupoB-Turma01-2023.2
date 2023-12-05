//=======================================================================================//
//                 PME3380 - Modelagem de Sistemas Dinâmicos - 2023.2                    //
//                  Análise de resposta no domínio do tempo - Pulso                      //
//=======================================================================================//
//                               Versão do Scilab: 5.5.2                                 //
//=======================================================================================//
//Eduardo Henrique de Azevedo Cabrita                                            12553462//
//Gabriel Circeli Barbieri Santin                                                12624530//
//Iago Cardoso Nogueira da Silva                                                 10738721//
//Lucas Sung Il Hong                                                             12717287//
//=======================================================================================//

//                                 SISTEMA NÃO LINEAR                                    //

//Constantes Numéricas:
clc;
clear;
L1=2.4;
L2=10.9;
g=9.81;
m1=7030;
m2=3500;
IG1=39888;
IG2=1242;
c=10000;
rho=1.225;
A1=11.9;
A2=5.6;
C1=1.24;
C2=0.47;

//Espaço de Estados não-Linear:
funcprot(0);
    function dx=EENL(t,x)
    //Vetor de entradas (u):
    t=linspace(0,100,10000);
    for i=1:10000
        if t(i)<=33.3 then
        u=[0;0];
        elseif t(i)>33.3 & t(i)<=66.6
        u=[-60000;20000];
        else
        u=[0;0];
        end
    end
    
    dx(1)=x(4);
    
    dx(2)=x(5);
    
    dx(3)=(2*sec(x(2))*(((m2*L1^2+IG1)*cos(x(2))-L1^2*m2*cos(x(1))*cos(x(1)-x(2)))*(L2*m2*cos(x(2))*(g*L2*m2*sin(x(2))-x(4)*(c+L1*L2*m2*sin(x(1)-x(2))*x(4))+c*x(5)+1/2*L2*rho*cos(x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5)))-(m2*L2^2+IG2)*(-L1*m2*sin(x(1))*x(4)^2-L2*m2*sin(x(2))*x(5)^2-u(1)+1/2*rho*(abs(x(3))*C1*A1*x(3)+abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5)))))-L1*m2*(L2^2*m2*cos(x(1)-x(2))*cos(x(2))-(m2*L2^2+IG2)*cos(x(1)))*(cos(x(2))*(L1*L2*m2*sin(x(1)-x(2))*x(5)^2-c*x(5)-u(2)+g*L1*m2*sin(x(1))+c*x(4)+1/2*L1*rho*cos(x(1))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5)))-L1*cos(x(1)-x(2))*(-L1*m2*sin(x(1))*x(4)^2-L2*m2*sin(x(2))*x(5)^2-u(1)+1/2*rho*(abs(x(3))*C1*A1*x(3)+abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5)))))))/(IG2*m2^2*L1^2+L2^2*m1*m2^2*L1^2+2*IG2*m1*m2*L1^2-IG2*m2^2*cos(2*x(1))*L1^2-L2^2*m1*m2^2*cos(2*(x(1)-x(2)))*L1^2+IG1*L2^2*m2^2+2*IG1*IG2*m1+2*IG1*IG2*m2+2*IG1*L2^2*m1*m2-IG1*L2^2*m2^2*cos(2*x(2)));
    
    dx(4)=-((8*L1*m1*m2^2*sin(x(1)-x(2))*x(5)^2*L2^3-L1*m1*m2*rho*cos(x(1)-3*x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L2^3+L1*m1*m2*rho*cos(x(1)+x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L2^3+4*L1^2*m1*m2^2*sin(2*(x(1)-x(2)))*x(4)^2*L2^2+4*g*L1*m1*m2^2*sin(x(1))*L2^2+4*g*L1*m1*m2^2*sin(x(1)-2*x(2))*L2^2-2*L1*m2^2*rho*cos(x(1))*abs(x(3))*C1*A1*x(3)*L2^2+2*L1*m2^2*rho*cos(x(1)-2*x(2))*abs(x(3))*C1*A1*x(3)*L2^2+2*L1*m1*m2*rho*cos(x(1))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(3)*L2^2-2*L1*m1*m2*rho*cos(x(1)-2*x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(3)*L2^2+4*c*m2^2*x(4)*L2^2+8*c*m1*m2*x(4)*L2^2-4*c*m2^2*cos(2*x(2))*x(4)*L2^2+L1^2*m1*m2*rho*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L2^2+L1^2*m1*m2*rho*cos(2*x(1))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L2^2-L1^2*m1*m2*rho*cos(2*(x(1)-x(2)))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L2^2-L1^2*m1*m2*rho*cos(2*x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L2^2-4*c*m2^2*x(5)*L2^2-8*c*m1*m2*x(5)*L2^2+4*c*m2^2*cos(2*x(2))*x(5)*L2^2+4*IG2*L1*m2^2*sin(x(1)-x(2))*x(5)^2*L2+8*IG2*L1*m1*m2*sin(x(1)-x(2))*x(5)^2*L2+4*IG2*L1*m2^2*sin(x(1)+x(2))*x(5)^2*L2+4*c*L1*m2^2*cos(x(1)-x(2))*x(4)*L2+8*c*L1*m1*m2*cos(x(1)-x(2))*x(4)*L2-4*c*L1*m2^2*cos(x(1)+x(2))*x(4)*L2-4*c*L1*m2^2*cos(x(1)-x(2))*x(5)*L2-8*c*L1*m1*m2*cos(x(1)-x(2))*x(5)*L2+4*c*L1*m2^2*cos(x(1)+x(2))*x(5)*L2+2*IG2*L1*m1*rho*cos(x(1)-x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L2+2*IG2*L1*m1*rho*cos(x(1)+x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L2+4*IG2*L1^2*m2^2*sin(2*x(1))*x(4)^2-4*L1*m2*(L2^2*m2*cos(x(1)-2*x(2))-(m2*L2^2+2*IG2)*cos(x(1)))*u(1)-4*(m2*(2*m1+m2)*L2^2-m2^2*cos(2*x(2))*L2^2+2*IG2*(m1+m2))*u(2)+8*g*IG2*L1*m2^2*sin(x(1))+8*g*IG2*L1*m1*m2*sin(x(1))-4*IG2*L1*m2*rho*cos(x(1))*abs(x(3))*C1*A1*x(3)+4*IG2*L1*m1*rho*cos(x(1))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(3)+8*c*IG2*m1*x(4)+8*c*IG2*m2*x(4)+2*IG2*L1^2*m1*rho*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)+2*IG2*L1^2*m1*rho*cos(2*x(1))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)-8*c*IG2*m1*x(5)-8*c*IG2*m2*x(5))/(4*(IG2*m2^2*L1^2+L2^2*m1*m2^2*L1^2+2*IG2*m1*m2*L1^2-IG2*m2^2*cos(2*x(1))*L1^2-L2^2*m1*m2^2*cos(2*(x(1)-x(2)))*L1^2+IG1*L2^2*m2^2+2*IG1*IG2*m1+2*IG1*IG2*m2+2*IG1*L2^2*m1*m2-IG1*L2^2*m2^2*cos(2*x(2)))));
    
    dx(5)=-((-8*L2*m1*m2^2*sin(x(1)-x(2))*x(4)^2*L1^3-L2*m1*m2*rho*cos(3*x(1)-x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L1^3+L2*m1*m2*rho*cos(x(1)+x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L1^3-4*L2^2*m1*m2^2*sin(2*(x(1)-x(2)))*x(5)^2*L1^2-4*g*L2*m1*m2^2*sin(2*x(1)-x(2))*L1^2+4*g*L2*m1*m2^2*sin(x(2))*L1^2+2*L2*m2^2*rho*cos(2*x(1)-x(2))*abs(x(3))*C1*A1*x(3)*L1^2-2*L2*m2^2*rho*cos(x(2))*abs(x(3))*C1*A1*x(3)*L1^2-2*L2*m1*m2*rho*cos(2*x(1)-x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(3)*L1^2+2*L2*m1*m2*rho*cos(x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(3)*L1^2-4*c*m2^2*x(4)*L1^2-8*c*m1*m2*x(4)*L1^2+4*c*m2^2*cos(2*x(1))*x(4)*L1^2+4*c*m2^2*x(5)*L1^2+8*c*m1*m2*x(5)*L1^2-4*c*m2^2*cos(2*x(1))*x(5)*L1^2+L2^2*m1*m2*rho*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L1^2-L2^2*m1*m2*rho*cos(2*x(1))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L1^2-L2^2*m1*m2*rho*cos(2*(x(1)-x(2)))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L1^2+L2^2*m1*m2*rho*cos(2*x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)*L1^2-4*IG1*L2*m2^2*sin(x(1)-x(2))*x(4)^2*L1-8*IG1*L2*m1*m2*sin(x(1)-x(2))*x(4)^2*L1+4*IG1*L2*m2^2*sin(x(1)+x(2))*x(4)^2*L1+8*L2*m2*u(2)*(m1*cos(x(1))*cos(x(2))+(m1+m2)*sin(x(1))*sin(x(2)))*L1-4*c*L2*m2^2*cos(x(1)-x(2))*x(4)*L1-8*c*L2*m1*m2*cos(x(1)-x(2))*x(4)*L1+4*c*L2*m2^2*cos(x(1)+x(2))*x(4)*L1+2*IG1*L2*m1*rho*cos(x(1)-x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L1+2*IG1*L2*m1*rho*cos(x(1)+x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(4)*L1+4*c*L2*m2^2*cos(x(1)-x(2))*x(5)*L1+8*c*L2*m1*m2*cos(x(1)-x(2))*x(5)*L1-4*c*L2*m2^2*cos(x(1)+x(2))*x(5)*L1+4*IG1*L2^2*m2^2*sin(2*x(2))*x(5)^2-4*L2*m2*(L1^2*m2*cos(2*x(1)-x(2))-(m2*L1^2+2*IG1)*cos(x(2)))*u(1)+8*g*IG1*L2*m2^2*sin(x(2))+8*g*IG1*L2*m1*m2*sin(x(2))-4*IG1*L2*m2*rho*cos(x(2))*abs(x(3))*C1*A1*x(3)+4*IG1*L2*m1*rho*cos(x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(3)-8*c*IG1*m1*x(4)-8*c*IG1*m2*x(4)+8*c*IG1*m1*x(5)+8*c*IG1*m2*x(5)+2*IG1*L2^2*m1*rho*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5)+2*IG1*L2^2*m1*rho*cos(2*x(2))*abs(x(3)+L1*cos(x(1))*x(4)+L2*cos(x(2))*x(5))*C2*A2*x(5))/(4*(IG2*m2^2*L1^2+L2^2*m1*m2^2*L1^2+2*IG2*m1*m2*L1^2-IG2*m2^2*cos(2*x(1))*L1^2-L2^2*m1*m2^2*cos(2*(x(1)-x(2)))*L1^2+IG1*L2^2*m2^2+2*IG1*IG2*m1+2*IG1*IG2*m2+2*IG1*L2^2*m1*m2-IG1*L2^2*m2^2*cos(2*x(2)))));
endfunction;

//Condição Inicial:
theta1_0=0.492422;
theta2_0=0.276743;
x1p_0=-280/3.6;
theta1p_0=0;
theta2p_0=0;

//Solução:
t=linspace(0,100,10000);
sol=ode([theta1_0;theta2_0;x1p_0;theta1p_0;theta2p_0],0,t,EENL);
theta1=sol(1,:);
theta2=sol(2,:);
x1p=sol(3,:);

//Gráficos:
f1=scf(1)
plot(t,theta1,"b")
xtitle('Theta1 em função do tempo')

f1=scf(2)
plot(t,theta2,"b")
xtitle('Theta2 em função do tempo')

f1=scf(3)
plot(t,x1p,"b")
xtitle('X1p em função do tempo')

//                                 SISTEMA     LINEAR                                    //

//Espaço de Estados Linear:
funcprot(0);
    function dx=EEL(t,x)
    //Vetor de entradas (u):
    t=linspace(0,100,10000);
    for i=1:10000
        if t(i)<=33.3 then
        u=[0;0];
        elseif t(i)>33.3 & t(i)<=66.6
        u=[-60000;20000];
        else
        u=[0;0];
        end
    end
    A=[0,0,0,1,0;0,0,0,0,1;
    -0.131251,4.82312,-0.195317,-0.145171,0.112308;
    -2.04901,2.0966,-0.00141516,-0.297085,0.307039;
    0.450941,-1.8068,0.0111849,0.0871245,-0.165803];
    B=[0,0;
    0,0;
    0.000137069,1.56851*10^-6;
    1.56851*10^-6,0.0000244866;
    -0.000012397,-5.38896*10^-6];
    dx=A*(x-[0.492422;0.276743;-280/3.6;0;0])+B*(u-[-64426.8;18332.9]);
endfunction;

//Condições Inicial:
//IDEM

//Vetor de entradas (u):
//IDEM

//Solução:
t=linspace(0,100,10000);
sol=ode([theta1_0;theta2_0;x1p_0;theta1p_0;theta2p_0],0,t,EEL);
theta1=sol(1,:);
theta2=sol(2,:);
x1p=sol(3,:);

//Gráficos:
f1=scf(1)
plot(t,theta1,"r")
xtitle('Theta1 em função do tempo')
h1=legend(['Sistema não-Linear','Sistema Linear'],["b","r"],4)

f1=scf(2)
plot(t,theta2,"r")
xtitle('Theta2 em função do tempo')
h1=legend(['Sistema não-Linear','Sistema Linear'],["b","r"],4)

f1=scf(3)
plot(t,x1p,"r")
xtitle('X1p em função do tempo')
h1=legend(['Sistema não-Linear','Sistema Linear'],["b","r"],5)
