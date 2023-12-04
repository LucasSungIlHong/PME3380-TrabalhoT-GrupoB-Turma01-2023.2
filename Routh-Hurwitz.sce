clc;
clear;

//=======================================================================================//
//                 PME3380 - Modelagem de Sistemas Dinâmicos - 2023.2                    //
//                Método de Routh-Hurwitz - Cálculo da Tabela de Routh                   //
//=======================================================================================//
//                               Versão do Scilab: 5.5.2                                 //
//=======================================================================================//
//Eduardo Henrique de Azevedo Cabrita                                            12553462//
//Gabriel Circeli Barbieri Santin                                                12624530//
//Iago Cardoso Nogueira da Silva                                                 10738721//
//Lucas Sung Il Hong                                                             12717287//
//=======================================================================================//

s=poly(0,'s');

//Critério de Routh-Hurwitz:
[r,num]=routh_t(s^5+0.658204*s^4+3.96726*s^3+1.25885*s^2+2.85069*s+0.433712);
disp(r)

disp('trocas de sinal na primeira coluna da tabela.',num,'Há')
