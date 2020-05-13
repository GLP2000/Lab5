#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "winbgi2.h"
#include "rk4.h"

void rhs_fun(double t, double* G, double* F);
void veuler(double t, double* G, double h, int n,void (*fun)(double, double*, double*), double* G1);
double energia(double omega, double alfa);
double g = 9.8; 
double l = 0.5;
double m = 2.0;

int main() 
{
	double s, omega0;
	printf("podaj kat alfa0 w stopni\n");
    scanf("%lf", &s);
	printf("podaj omega0\n");
	scanf("%lf", &omega0);

	double  alfa0= s/180.*3.141593;
	
	double t;
	double t1 = 0;
	double t2 = 10;
	double h = 0.001;
	double G[2] = { omega0,alfa0 }, G1[2];
	t = t1;
	double ec;

	graphics(600, 400);
	scale(-5, -4, 5,4);
	title("Metoda Eulera(zielony) Metoda RK4(niebielski)", " ", " ");
    


	printf("\nLiczymy dla metody Eulery\n");

	while (t < t2)
	{

		veuler(t, G, h, 2, rhs_fun, G1);

		G[0] = G1[0];
		G[1] = G1[1];

		t += h;

		ec = energia(G[0], G[1]);
		printf("Omega=%lf\tAlfa=%lf\tEnergy=%lf\n", G[0], G[1], ec);

		setcolor(0.5);//zielone
		point(G[0], G[1] );

		

	}
	t = t1;
	G[0] = omega0;
	G[1] = alfa0;

	printf("\nLiczymy dla metody RK4\n");

	while ( t < t2)
	{
		vrk4(t, G, h, 2, rhs_fun, G1);

		G[0] = G1[0];
		G[1] = G1[1]; 

		t+= h;

		setcolor(0.125);//niebielski
		point(G[0], G[1]);

		ec = energia(G[0], G[1]);
		printf("Omega RK4=%lf\tAlfa RK4=%lf\tEnergy RK4=%lf\n", G[0], G[1], ec);

		
	}



	wait();


}

void rhs_fun(double t, double* G, double* F)
{

	F[1] = -(g / l) * sin(G[0]);
	F[0] = G[1];

}

void veuler(double t, double* G, double h, int n, void (*fun)(double, double*, double*), double* G1)
{
	fun(t, G, G1);
	for (int i = 0; i < n; i++) 
	{
		G1[i] = G[i] + h * G1[i];

	}

}

double energia(double omega, double alfa)
{
	return (m * l * l) / 2 * (omega*omega) + m * g * l * (1 - cos(alfa));
}
