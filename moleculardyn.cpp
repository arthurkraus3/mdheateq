#include <stdio.h>
#include <iostream>
#include <math.h>

#include "fileutils.h"
#include "stringutils.h"
#include "gnuplot.h"

#define g0 9.8067 //default gravity

#define eps 1e-8

//LJ parameters
#define epsilon 1.0
#define sigma 1.0


typedef struct {
	double x,y;
} point2D;

typedef struct {
    double T; //system temperature
    double g; //gravity
    int BC; //bit flag: 1: peridoc in x-direction if 1, bit 2: open at top if 1, i.e., possbile values: 0 (closed),1 (periodic in x), 2 (open at top), 3 (open and peridic) --- open at top only used if g>0
} sysparams;

//---------------------------------------------------------------------------

//calculate LF forces between particle i and j using const LF parameters above, distances are cut-off by eps
void LJforce(point2D *r,int i,int j,point2D &f) {
    point2D dr;
    double r2,ir2,eir2,a,x;
    //printf("Initialized calcForces()");

    //i & j are the particle indices, f.x and f.y the x and y components of the force. Here we just exclude self-interation and return with a zero force
    if(i==j) {f.x=f.y=0.0;return;}
    
    //distance calculation in x and y direction, note the sign!
    dr.x=r[i].x-r[j].x;
    dr.y=r[i].y-r[j].y;
    
    //calculate the eucledean distance (squared), if too small cut-off at eps (here 1e-8)
    r2=dr.x*dr.x+dr.y*dr.y;
    if(r2<eps) r2=eps;
    
    //inverse distance and force calculation
    ir2=1.0/r2;
    eir2=epsilon*epsilon*ir2;
    a=eir2*eir2*eir2;
    x=24.0*sigma*ir2*a*(2.0*a-1.0);
    //again dr.x and dr.y determine the sign of the force
    f.x=x*dr.x;
    f.y=x*dr.y;
}

//---------------------------------------------------------------------------

//calculate all forces
void calcforces(point2D *r,point2D *f,double m,int N,double g=g0) {
    int i,j;
    point2D tf;

    //printf("Initialized calcForces()");
    
    for(i=0;i<N;i++) {
        f[i].x=0.0;
        f[i].y=-m*g;
    } //initialize all forces by gravitational force
    
    for(i=0;i<N;i++) { //add all pair forces
        for(j=i+1;j<N;j++) {
            LJforce(r,i,j,tf);
            f[i].x+=tf.x;
            f[i].y+=tf.y;
            f[j].x-=tf.x;
            f[j].y-=tf.y;
        }}
    
}
//---------------------------------------------------------------------------

inline void reflectbox(point2D &r,point2D &v,double Lx,double Ly,int BC=0) {
    //x direction
    if(r.x>Lx) {
        if((BC&1)==1) {r.x-=Lx;} //periodic
        else {r.x=Lx-(r.x-Lx);v.x=-v.x;} //reflection at the right boundary
    }
    else if(r.x<0) {
        if((BC&1)==1) {r.x+=Lx;} //periodic
        else {r.x=-r.x;v.x=-v.x;} //reflection at the left boundary
    }
    
    //y direction
    if((r.y>Ly) && ((BC&2)==0)) {
        r.y=Ly-(r.y-Ly);v.y=-v.y;
    }
    else if(r.y<0) {
        r.y=-r.y;
        v.y=-v.y;    
    } //always reflection at bottom (if gravity is zero one could also make it periodic)
}

//---------------------------------------------------------------------------

//velocity verlet
void velverlet(point2D *r,point2D *v,point2D *f,point2D *f1,double m,double ht,double Lx,double Ly,int N,sysparams &sc) {
    int i;
    point2D org;
    double T,v2;
    //printf("Initialized velverlet()");
    //update positions
    for(i=0;i<N;i++) {
        org=r[i];
        r[i].x=org.x+(v[i].x+0.5*f[i].x*ht)*ht;
        r[i].y=org.y+(v[i].y+0.5*f[i].y*ht)*ht;
        reflectbox(r[i],v[i],Lx,Ly,sc.BC);
    }
    calcforces(r,f1,m,N,sc.g); //calculate forces for new positions - reuse for next timestep
    
    T=0.0;
    //update velocities and calculate temperature (also in the left and right halves (even without wall))
    for(i=0;i<N;i++) {
        v[i].x=v[i].x+0.5*(f[i].x+f1[i].x)*ht;
        v[i].y=v[i].y+0.5*(f[i].y+f1[i].y)*ht;
        v2=v[i].x*v[i].x+v[i].y*v[i].y;
        T+=v2;
    }

    sc.T=m*T/(2.0*(N-1));
}

//---------------------------------------------------------------------------
//outputs a ascii file with x,y,vx,vy per particle
//creates gnuplot script file and calls gnuplot with this script to generate png file

void output(point2D *r,point2D *v,double m,double Lx,double Ly,int N,double t,sysparams &sc,int fn,gnuplot *gp=NULL) {
    string s,fname;
    fHandle f;
    int i;
    double a;
    //printf("Initialized output()\n");
    
    //output
    fname="OUT/MD_"+IntToStrF(fn,6)+".txt";
    f=FileCreate(fname);
    //printf("output(): created f\n");

    for(i=0;i<N;i++) {
        //printf("output(): for %i \n",i);
        s=FloatToStr(r[i].x)+"\t"+FloatToStr(r[i].y)+"\t"+FloatToStr(v[i].x)+"\t"+FloatToStr(v[i].y)+"\n";
        //std::cout << s << "\n";
        FileWrite(f,s.c_str(),s.length());
    }
    //printf("output(): File Written\n");
    
    FileClose(f);
    
    if(gp==NULL) return;
    
    //plot
    a=Ly/Lx; //aspect ratio
    gp->setterm("png medium size 1280,"+IntToStr((int) (1280*a))+" enhanced");
    gp->setout("OUT/MD_"+IntToStrF(fn,6)+".png");
    gp->addcommand("set grid xtics\nset grid ytics\nshow grid");
    gp->addcommand("set title 'time: "+FloatToStr(t)+", T="+FloatToStr(sc.T)+"'");
    gp->plotfile(fname,"u 1:2:(0.2) w circles fill solid t ''","[0:"+FloatToStr(Lx)+"] [0:"+FloatToStr(Ly)+"]");
    gp->setout();
    gp->run("2> /dev/null");
    gp->reset();
    
    printf("%d\t%le\t%le\n",fn,t,sc.T);
}


//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


int main() {
    //initial variables
	int N=100;
    int Nt=5000;
    double Lx=30;
    double Ly=30;
    double m=1.0;
    double ht=0.001;
    double v0=5.0; //initial max velocity
    
    int n,i,fn,nc;
    double t;
    double x,y;
    point2D *r,*v,*f,*f1,*ft;
    sysparams scc,sca;
    
    scc.g=0.0; //g0;
    scc.BC=2; //open top: 2

    //printf("Initialized main()\n");
    
    //avoid open top without gravity
    if((scc.g<eps) && ((scc.BC&2)==2)) scc.BC&=1;//&= makes scc.BC = 1
    
    
    gnuplot *gp=new gnuplot();

    r=new point2D[N]; //position
    v=new point2D[N]; //velocity
    f=new point2D[N]; 
    f1=new point2D[N];

    scc.T=0.0; //start with temperature = 0
    for(i=0;i<N;i++) { //iterate and assign values for r[i] and v[i] x and y components
        r[i].x=10.5+1.0*(i%10);
        r[i].y=10+1.0*(i/10);
        v[i].x=v0*(2.0*drand48()-1.0);
        v[i].y=v0*(2.0*drand48()-1.0);
        scc.T+=v[i].x*v[i].x+v[i].y*v[i].y;
    }
    scc.T=m*scc.T/(2.0*(N-1));//T = 1/2 m vx^2 + vy^2

    fn=0;
    
    output(r,v,m,Lx,Ly,N,0.0,scc,fn,gp);fn++;
    
    sca=scc; //sca used for averaging of temperatures and particle numbers
    nc=1;
    calcforces(r,f,m,N,scc.g); //calculate initial forces
	for(n=0;n<Nt;n++) {
        t=n*ht;
        velverlet(r,v,f,f1,m,ht,Lx,Ly,N,scc);
        
        //swap forces
        ft=f1;f1=f;f=ft;
        
        sca.T+=scc.T;

        nc++;
        
        if((n+1)%100==0) {
            sca.T=sca.T/nc;
            output(r,v,m,Lx,Ly,N,t+ht,sca,fn,gp);fn++;
            nc=0;
            sca.T=0.0;
        }
	}
    //create movie from png images created by gnuplot
    system("ffmpeg -y -i OUT/MD_%06d.png OUT/MD.m4v");
	
	delete[] r;
	delete[] v;
    delete[] f;
    delete[] f1;
	delete gp;
	return 0;
}
