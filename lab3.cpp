#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include "math.h"
using namespace std;
ofstream archivo;
#define pi 3.141592654
#define ep 2.718281828

double fun( double x){
	return (1/(pow(3,x)));
	}
double fun2(double x){
	return (2/pow(3,pow(2,x)));
}
double fun3(double x){
	return (1/pow(x,2));
}
double fun4(double x){
	if(x>=0){
		if(x==0)
			return 2;
		else{
			return ((1/2*fun4(x-1))+(1/fun4(x-1)));
		}
}
}
double fun5(double x){
	return pow(x,3)+4*pow(x,2)-10;
}
double d_fun5(double x){
	float h=0.0001;
	return ((fun5(x+h)-fun(x))/h);
}
double newton(double r,double a,double b){
	if(r>=a && r<=b)
		return r-(fun5(r)/d_fun5(r));
}
double alpha_1(int n){
	return ((log(abs(fun(n-1)))-log(abs(fun(n))))/(log(abs(fun(n)))/log(abs(fun(n-1)))));
}
double alpha_2(int n,double x0=0){
	return ((log(abs(fun2(n-1)-x0))-log(abs(fun2(n)-x0)))/(log(abs(fun2(n)-x0))/log(abs(fun2(n-1)-x0))));
}
double alpha_3(int n,double x0=0){
	return ((log(abs(fun3(n-1)-x0))-log(abs(fun3(n)-x0)))/(log(abs(fun3(n)-x0))/log(abs(fun3(n-1)-x0))));
}
double alpha_4(int n,double x0=0){
	return ((log(abs(fun4(n-1)-x0))-log(abs(fun4(n)-x0)))/(log(abs(fun4(n)-x0))/log(abs(fun4(n-1)-x0))));
}


/*double falsaposicion_tolerancia(double funcion(double X), double parametro1, double parametro2, double tolerancia){
	archivo.open("datos.txt");
	archivo<< "i" << "\t" << "parametro1" << "\t" << "parametro2" << " \t" << "iteracion"<< " \t" << "funcion"<<endl;
	double m = parametro2-funcion(parametro2)*((parametro2-parametro1)/(funcion(parametro2)-funcion(parametro1)));
	int i=0;
	while(funcion(m)>tolerancia || funcion(m)<-tolerancia){
		archivo<< i << " " << parametro1 << " " << parametro2 << " " << tolerancia << " " << funcion(tolerancia)<<endl;
		if (funcion(parametro1)*funcion(m)<0){
			parametro2=m;
		}
		else{
			parametro1=m;
		}
		i++;
		m = parametro2-funcion(parametro2)*((parametro2-parametro1)/(funcion(parametro2)-funcion(parametro1)));

	}
	return m;
}*/

int main(){
	double i=1;
	int n=2;
	double r0;
	double r=1.5;
 	while(n<100){
 		/*i=abs(alpha_4(n)-alpha_4(n-1));
 		if(i>0.0001){
	 		cout<<n-1<<":"<<alpha_4(n)<<" Error: "<<i<<endl;
	 		n++;
	 	}	else break;*/
	 	//NEWTON
	 	r0=r;
	 	n++;
	 	r=newton(r,1,2);
	 	i=abs(r-r0);
	 	cout<<n-1<<": "<<r<<" Error: "<<i<<endl;



 	}
	//archivo.close();
	return 0;
	}
