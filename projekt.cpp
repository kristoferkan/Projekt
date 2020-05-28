#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "rk4.h"
#include "winbgi2.h"
double h;
double x0;
double v0;
double m;
double k1;
double k2;
double c;

void rhs_fun(double t,double *X, double *F);
void energia(double *X,double *Ec, double *Ep, double *Ek);

void main()
{
	double Ec,Ep,Ek;
	double X[2];
	double F[2];
	double X1[2];
	int n;
	printf("Podaj krok calkowania:\n");
	scanf("%lf",&h);
	printf("Podaj poczatkowe wychylenie:\n");
	scanf("%lf",&x0);
	printf("Podaj poczatkowa predkosc:\n");
	scanf("%lf",&v0);
	printf("Podaj mase:\n");
	scanf("%lf",&m);
	printf("Podaj wspolczynnik sprezystosci k1:\n");
	scanf("%lf",&k1);
	printf("Podaj wspolczynnik sprezystosci k2:\n");
	scanf("%lf",&k2);
	printf("Podaj wspolczynnik tlumienia:\n");
	scanf("%lf",&c);
	printf("Jaki wykres wyswietlic:\nx(t)-0\nv(t)-1\nE(t)-2\nfazowy-3\n");
	scanf("%d",&n);
	graphics(800, 600);
	X[0]=x0;
	X[1]=v0;
	switch(n)
	{
		case 0:
			scale(0, -2, 20, 2);
			title("t","x","");
			for(double t=0.;t<=20.;t+=h)
			{	
				setgray(1.);
				point(t,0);
				setcolor(0.);
				point(t,X[0]);
				vrk4(t,X,h,2,rhs_fun,X1);
				X[0]=X1[0];
				X[1]=X1[1];
			}
			break;
		case 1:
			scale(0, -6, 20, 6);
			title("t","v","");
			for(double t=0.;t<=20.;t+=h)
			{
				setgray(1.);
				point(t,0);
				setcolor(0.);
				point(t,X[1]);
				vrk4(t,X,h,2,rhs_fun,X1);
				X[0]=X1[0];
				X[1]=X1[1];
			}
			break;
		case 2:
			scale(0, 0, 20, 20);
			title("t","E","Energia calkowita ukladu");
			for(double t=0.;t<=20.;t+=h)
			{
				energia(X,&Ec,&Ep,&Ek);
				setcolor(0.5);
				point(t,Ec);
				vrk4(t,X,h,2,rhs_fun,X1);
				X[0]=X1[0];
				X[1]=X1[1];
			}
			break;
		case 3:
			scale(-1.5, -6, 1.5, 6);
			title("x","v","");
			for(double t=0.;t<=20.;t+=h)
			{
				setcolor(0.);
				point(X[0],X[1]);
				vrk4(t,X,h,2,rhs_fun,X1);
				X[0]=X1[0];
				X[1]=X1[1];
			}
			break;
	}

	wait();
}
void rhs_fun(double t,double *X, double *F)
{
	F[0] = X[1];
	F[1] =-1/m*(k1*(1+k2*pow(X[0],2))*X[0]+c*X[1]);
}
void energia(double *X,double *Ec, double *Ep, double *Ek)
{
	*Ep=m*pow(X[1],2)/2;
	*Ek=k1*pow(X[0],2)/2+k2*pow(X[0],4)/4;
	*Ec=*Ep+*Ek;
}