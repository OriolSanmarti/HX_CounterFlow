clear all
close all

[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp==6.3.0'])
fluid='Water';

%Dades entrada
T1i = 100;
T2i = 25;
Di = 0.014;
Do = 0.014;
De = 0.025;
S1 = (Di/2)^2 * pi;
S2 = (( (De^2) - (Do^2) ) / 4) * pi;
alphai = 0;
alphao = 0;
Rfi = 0;%5E-5;%5E-3;
Rfo = 0;%5E-5;%5E-3;
lambda = 430;
L = 10;
p1i = 6.0e6;
p2i = 6.0e6;
error = 0.0001;
%ro1  = py.CoolProp.CoolProp.PropsSI('D', 'T', T1i+273.15, 'P', p1i, fluid);
%ro2  = py.CoolProp.CoolProp.PropsSI('D', 'T', T2i+273.15, 'P', p2i, fluid);
ro1 = 1000;
ro2 = 1000;
m1 = 0.46;
m2 = 0.8466;
v1i = m1 / (S1 * ro1);
v2i = m2 / (S2 * ro2);

%Step 0 Estimate properties
T1o = 40;
T2o = 50;
p1o = 6.0e6;
p2o = 6.0e6;
Tw1 = 85;
Tw2 = 80;

itefin = false;
while itefin == false;
    itefin = true;
    %Step 1 evaluate thermo-properties
    T1 = 273.15 + ((T1i + T1o) / 2);
    T2 = 273.15 + ((T2i + T2o) / 2);
    p1 = (p1i + p1o) / 2;
    p2 = (p2i + p2o) / 2;
    Tw1ref = Tw1 + 273.15;
    Tw2ref = Tw2 + 273.15;
    ro1  = py.CoolProp.CoolProp.PropsSI('D', 'T', T1, 'P', p1, fluid);
    cp1  = py.CoolProp.CoolProp.PropsSI('C', 'T', T1, 'P', p1, fluid);
    mu1  = py.CoolProp.CoolProp.PropsSI('V', 'T', T1, 'P', p1, fluid);
    muw1  = py.CoolProp.CoolProp.PropsSI('V', 'T', Tw1ref, 'P', p1, fluid);
    ro2  = py.CoolProp.CoolProp.PropsSI('D', 'T', T2, 'P', p2, fluid);
    cp2  = py.CoolProp.CoolProp.PropsSI('C', 'T', T2, 'P', p2, fluid);
    mu2  = py.CoolProp.CoolProp.PropsSI('V', 'T', T2, 'P', p2, fluid);
    muw2  = py.CoolProp.CoolProp.PropsSI('V', 'T', Tw2ref, 'P', p2, fluid);
    lambda1 = py.CoolProp.CoolProp.PropsSI('L', 'T', T1, 'P', p1, fluid);
    lambda2 = py.CoolProp.CoolProp.PropsSI('L', 'T', T2, 'P', p2, fluid);
    v1o = m1/(ro1*S1);
    v2o = m2/(ro2*S2);
    v1 = (v1o + v1i) / 2;
    v2 = (v2i + v2o) / 2;
    Re1 = (ro1 * v1 * Di) / (mu1);
    Re2 = (ro2 * v2 * (De-Do)) / (mu2);
    
    %Step 2 Calculate Uo
        %Calcul coeficcients conveccio
        
        Pr1 = (cp1 * mu1) / lambda1;
        Nu1 = 0.027 * Re1^0.8 * Pr1^0.33 * (mu1/muw1)^0.14;
        alphai = Nu1*lambda1/Di;
        Pr2 = (cp2 * mu2) / lambda2;
        Nu2 = 0.027 * Re2^0.8 * Pr2^0.33 * (mu2/muw2)^0.14;
        alphao = Nu2*lambda2/(De-Do);
    
    Po = pi * Do;
    Pi = pi * Di;
    Uo = ( ( ((1/alphai) + (Rfi)) * (Po / Pi) ) + ( (Po / (2 * pi * lambda)) * (log(Do/Di)) ) + (Rfo + (1/alphao)) )^(-1);
    Ao = L * (pi * Do);
    Ai = L * (pi * Di);
    
    %Step 3 Calculate epsilon
    C1 = m1 * cp1;
    C2 = m2 * cp2;
    Cmin = C1;     
    Cmax = C2;
    if C2 < C1
        Cmin = C2;
        Cmax = C1;
    end
    Z = Cmin / Cmax;
    NTU = (Uo * Ao) / (Cmin);
    %epsilon = (1 - exp( -NTU * (1 + Z)) ) / (1 + Z);
    epsilon = (1 - exp( -NTU * (1 - Z)) ) / (1 - Z*exp( -NTU * (1 - Z)) );
    
    %Step 4 Calculate Q
    Q = epsilon * Cmin * (T1i - T2i);

    %Step 5 Calculate T1o T2o
    T1oant = T1o;
    T1o = T1i - (Q / (m1 * cp1));
    T2oant = T2o;
    T2o = T2i + (Q / (m2 * cp2)); 
    if abs(T1o-T1oant) > error || abs(T2o-T2oant) > error
        itefin = false;
    end
    
    %Step 6 Evaluate Tw1 Tw2
    Tw1ant = Tw1;
    T1 = (T1o + T1i) / 2;
    Tw1 = T1 - ( (Q * (Rfi + (1/alphai))) / (Ai) );
    
    Tw2ant = Tw2;
    T2 = (T2o + T2i) / 2;
    Tw2 = T2 + ( (Q * (Rfo + (1/alphao))) / (Ao) );
    
    %Step 7 Calculate p1 p2
    rug = 0.0013/1000;
    
    er1 = rug/Di;
    f1 = 0.0625 / (log10(er1/3.7 + 5.74/(Re1^0.9)))^2;
    tau1 = f1*ro1*v1^2 / 2;
    incp1 = ( (tau1*pi*Di*L) + (m1*(v1o-v1i)) ) / (pi*Di^2/4);
    p1o = p1i - incp1;

    er2 = rug/De;
    f2 = 0.0625 / (log10(er2/3.7 + 5.74/(Re2^0.9)))^2;
    tau2 = f2*ro2*v2^2 / 2;
    incp2 = ((tau2*pi*(De+Do)*L) + (m2*(v2o-v2i)) ) / (pi*(De^2-Do^2)/4);
    p2o = p2i - incp2;
    
end
incp1 = p1o-p1i;
incp2 = p2o-p2i;
Pot1 = m1*cp1*(T1o-T1i)
Pot2 = m2*cp2*(T2o-T2i)
