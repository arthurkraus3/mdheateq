#include <stdio.h>
#include <iostream>
#include <math.h>

#include "fileutils.h"
#include "stringutils.h"
#include "gnuplot.h"


double Gamma(double Theta, double l, double L, double x){
    double il,y;

    il = 1.0/l;
    y = x-0.5*L;

    return (Theta*il*exp(-y*y*il*il));
}


int main() {
    int N = 1000;
    double L = 10.0;//larger L, more dropoff
    double l = 1.0;
    double kappa = 1.0;
    double Theta = -0.4;
    double x,g,h,beta,T1;
    double *T,*alpha;
    int i;

    string s;
    fHandle f;
    gnuplot *gp = new gnuplot();

    T = new double[N+1];
    alpha = new double[N+1];

    T[0] = 0.0;//fixed
    T[N] = 2.0;//fixed
    h = L / (1.0*N);

    //Solve the heat equation

    alpha[1] = -2.0; //from second derivative--
    beta = 0.0;
    for(i = 1; i < N; i++){
        x = i*h;
        g=Gamma(Theta, l, L, x);
        if(i == 1) g -= kappa*T[0]/(h*h);
        else if(i==N-1) g -= kappa*T[N]/(h*h);
        g=g*h*h/kappa;
        // T[i] = T[0] + (T[N] - T[0])*i*h/L;//first step from slides
        T[i] = g-beta*T[i-1];//this is our y solution, starting i=1 this allows us to get a nonzero T[0] value
        beta=1.0/alpha[i];
        alpha[i+1] = -2.0 - beta;

        //now we need to solve for our Ux = y
    }
    // algo for previous loop//thomas algorithm Ly = f
    // T[1] = g;
    // T[2] = g - beta*T[1];
    // T[3] = g - beta*T[2]

    T1 = 0.0;
    for(i = N-1; i >= 1; i--){
        T1 = (T[i] - T1)/alpha[i];
        T[i] = T1;
    }
    //algo for previous loop
    // T[N-1] = T[N-1]/alpha[N-1]
    // T[N-2]=(T[N-2] - T[N-1])/alpha[N-2]
    // T[N-3]=(T[N-3] - T[N-2])/alpha[N-3]

    //output
    f = FileCreate("inhom_heat.txt");
    for(i = 0; i<N; i++){
        x = i*h;
        g = Gamma(Theta,l, L , x);
        s=FloatToStr(x)+"\t"+FloatToStr(T[i])+ "\t"+FloatToStr(g)+"\n";
        FileWrite(f,s.c_str(),s.length());
    }
    FileClose(f);
    

    //plot
    gp->setxylable("'x'","'T'");
    gp->plotfile("inhom_heat.txt","u 1:2 w l t 'T'");
    gp->replotfile("inhom_heat.txt","u 1:3 w l t 'G'");
    gp->show();
 
    delete[] T;
    delete[] alpha;
    delete gp;
    return 0;
}
