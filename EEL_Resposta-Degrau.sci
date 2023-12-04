//=======================================================================================//
//                 PME3380 - Modelagem de Sistemas Dinâmicos - 2023.2                    //
//                  Análise de resposta no domínio do tempo - Degrau                     //
//=======================================================================================//
//                               Versão do Scilab: 5.5.2                                 //
//=======================================================================================//
//Eduardo Henrique de Azevedo Cabrita                                            12553462//
//Gabriel Circeli Barbieri Santin                                                12624530//
//Iago Cardoso Nogueira da Silva                                                 10738721//
//Lucas Sung Il Hong                                                             12717287//
//=======================================================================================//

//Espaço de Estados Linear:
t=linspace(0,100,10000);
funcprot(0);
    function dx=EEL(t,x)
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

//Condição Inicial:
theta1_0=0.492422;
theta2_0=0.276743;
x1p_0=-280/3.6;
theta1p_0=0;
theta2p_0=0;

//Vetor de entradas (u):
u=[-52000;20000];

//Solução:
sol=ode([theta1_0;theta2_0;x1p_0;theta1p_0;theta2p_0],0,t,EEL);
theta1=sol(1,:);
theta2=sol(2,:);
x1p=sol(3,:);

//Gráficos:
f1=scf(1)
plot(t,theta1,"b")
xtitle('Theta1 em função do tempo')
h1=legend(['Sistema Linear'])

f1=scf(2)
plot(t,theta2,"r")
xtitle('Theta2 em função do tempo')
h1=legend(['Sistema Linear'])

f1=scf(3)
plot(t,x1p,"g")
xtitle('X1p em função do tempo')
h1=legend(['Sistema Linear'])
