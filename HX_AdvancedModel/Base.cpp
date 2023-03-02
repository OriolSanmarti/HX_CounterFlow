#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double CalcRho(double T){
    return 847.2 + 1.298*T - 2.657E-3*pow(T,2);
}
double CalcCp(double T){
    return 5648.8 - 9.140*T + 14.1E-3*pow(T,2);
}
double CalcMu(double T){
    return exp(7.867 - 0.077*T + 9.04E-5*T*T);
}
double CalcLambda(double T){
    return -0.722 + 7.168E-3*T - 9.137E-6*pow(T,2);
}

int main(){
    int N = 1000;
    double T1ini = 100;
    double T2ini = 25;
    double P1ini = 6E6;
    double P2ini = 6E6;
    vector <double> T1(N, T1ini);
    vector <double> T2(N, T2ini);
    vector <double> T1ant(N, T1ini);
    vector <double> T2ant(N, T2ini);
    vector <double> P1(N, P1ini);
    vector <double> P2(N, P2ini);
    vector <double> v1(N, 0.0);
    vector <double> v2(N, 0.0);
    double T1in = 100;
    double T2in = 25;
    double Tw1 = T1ini;
    double Tw2 = T2ini;

    //Inlet properties
    double m1 = 0.46;

    //Outlet properties
    double m2 = -0.8466;

    //Geometric properties
    double Di = 0.014;
    double Do = 0.014;
    double De = 0.025;

    double Pi = Di * M_PI;
    double Po = Do * M_PI;
    double Pe = De * M_PI;

    double S1 = M_PI * Di * Di / 4.0;
    double S2 = M_PI * ((De * De) - (Do * Do)) / 4.0;
    cout << "S1: " << S1 << " S2: " << S2 << endl;

    double L = 10;
    double incx = L/N;
    double UO = 0;
    double Ai = incx * Pi;
    double Ao = incx * Po;
    cout << "Ai: " << Ai << " Ao: " << Ao << endl;

    //Initialize properties
    double rho1 = CalcRho(T1ini+273.15);
    double cp1 = CalcCp(T1ini+273.15);
    double mu1 = CalcMu(T1ini+273.15);
    double muw1 = CalcMu(T1ini+273.15);
    double lambda1 = CalcLambda(T1ini+273.15);
    double rho2 = CalcRho(T2ini+273.15);
    double cp2 = CalcCp(T2ini+273.15);
    double mu2 = CalcMu(T2ini+273.15);
    double muw2 = CalcMu(T2ini+273.15);
    double lambda2 = CalcLambda(T2ini+273.15);

    double lambdas = 40;

    double Pr1 = 0, Pr2 = 0, Nu1 = 0, Nu2 = 0;
    double alphai = 0, alphao = 0, Re1 = 0, Re2 = 0;

    //Fouling factors
    double Rfi = 0;
    double Rfo = 0;
    double rug = 5E-5;

    cout << endl;

    double error = 1E-8;
    bool itefin = false;
    int ite = 0;
    while(!itefin){
    itefin = true;
        for(int i = 0; i<N; i++){
            double Te = 0;
            double Tw = 0;

            //Calc variable properties
            rho1 = CalcRho(T1[i]+273.15);
            cp1 = CalcCp(T1[i]+273.15);
            mu1 = CalcMu(T1[i]+273.15);
            muw1 = CalcMu(Tw1+273.15);
            lambda1 = CalcLambda(T1[i]+273.15);
            rho2 = CalcRho(T2[i]+273.15);
            cp2 = CalcCp(T2[i]+273.15);
            mu2 = CalcMu(T2[i]+273.15);
            muw2 = CalcMu(Tw2+273.15);
            lambda2 = CalcLambda(T2[i]+273.15);

            //Calc Convective term
            v1[i] = m1/(rho1*S1);
            v2[i] = m2/(rho2*S2);

            Re1 = (rho1 * v1[i] * Di) / mu1;
            Re2 = -(rho2 * v2[i] * (De - Do)) / mu2;

            Pr1 = (cp1 * mu1) / lambda1;
            Nu1 = 0.027 * pow(Re1, 0.8) * pow(Pr1, 0.33) * pow((mu1/muw1),0.14);
            alphai = Nu1 * lambda1 / Di;

            Pr2 = (cp2 * mu2) / lambda2;
            Nu2 = 0.027 * pow(Re2, 0.8) * pow(Pr2, 0.33) * pow((mu2/muw2),0.14);
            alphao = Nu2 * lambda2 / (De - Do);
            
            //Calculate UO
            UO = pow( ((((1/alphai) + Rfi) * (Po/Pi) + ( (Po / (2*M_PI*lambdas)) * (log(Do/Di)) ) + (Rfo + (1/alphao)))  ), -1);

            //Inlet
            if(i == 0) Tw = T1in;
            else Tw = T1[i-1];
            double Afi = m1*cp1 + UO*Po*incx;
            double Bfi = -m1*cp1;
            double Cfi = UO * Po * incx * T2[i];

            T1[i] = (Cfi - Bfi*Tw) / Afi;

            //Outlet
            if(i == N-1) Te = T2in;
            else Te = T2[i+1];
            double Afo = -m2*cp2 + UO*Po*incx;
            double Bfo = m2*cp2;
            double Cfo = UO * Po * incx * T1[i];

            T2[i] = (Cfo - Bfo*Te) / Afo;

            //Calculate Twall
            double Q = UO * (T1[i] - T2[i]) * Po * incx;
            Tw1 = T1[i] - ((Q * (Rfi + (1/alphai))) / (Ai));
            Tw2 = T2[i] + ((Q * (Rfo + (1/alphao))) / (Ao));

            //Calculate Pressure lost
                //Inlet
            double er1 = rug/Di;
            double ff1 = 0.0625 / pow(log10(er1/3.7 + 5.74/(pow(Re1, 0.9))), 2.0);
            double tau1 = ff1 * rho1 * v1[i] * v1[i] / 2.0;
            if(i == 0) P1[i] = P1ini;
            else P1[i] = P1[i-1] + (m1*(v1[i]-v1[i-1]) + tau1*Pi*incx) / S1;
                //Outlet
            double er2 = rug/De;
            double ff2 = 0.0625 / pow(log10(er2/3.7 + 5.74/(pow(Re2, 0.9))), 2.0);
            double tau2 = ff2 * rho2 * v2[i] * v2[i] / 2.0;
            if(i == N-1) P2[i] = P2ini;
            else P2[i] = P2[i+1] + (m2*(v2[i]-v2[i+1]) + tau2*Pe*incx) / S2;

            //Check error
            if(fabs(T2ant[i] - T2[i]) > error) itefin = false;
            if(fabs(T1ant[i] - T1[i]) > error) itefin = false;
            T1ant[i] = T1[i];
            T2ant[i] = T2[i];
        }
        
        ite++;

        cout << "\r" << "ite: " << ite<< " T1o: " << T1[N-1] << " T2o: " << T2[0] << " P1o: " << P1[N-1] << " P2o: " << P2[0]<< " incP1: " << P1[N-1] - P1[0] << " incP2: " << P2[0] - P2[N-1] << flush;

    }
    double Pot1 = 0.0, Pot2 = 0.0;
    for(int i = 1; i<T1.size(); i++){
        Pot1 += m1 * CalcCp(T1[i]+273.15) * (T1[i]-T1[i-1]);
        Pot2 += m2 * CalcCp(T2[i]+273.15) * (T2[i]-T2[i-1]);
    }
    cout << endl;
    cout << "Pot1: " << Pot1 << " Pot2: " << Pot2 <<endl;
    cout << endl << endl;
}
