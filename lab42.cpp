//https://juncotic.com/eliminacion-gaussiana-algoritmos-antiguos/
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include "math.h"
using namespace std;
ofstream archivo;
#define pi 3.141592654
#define ep 2.718281828
#define n 4

/*int** escalonado(int **A,int n){
}*/
void imprimir(double A[n][n+1]){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
		cout<<endl;
	}
}
/*void llenar(double &A[n][n+1]){
    int i,a;
	for(i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<i<<";"<<j<<": ";
			cin>>a;
			A[i][j]=a;
		}
		cout<<endl;
	}
}*/

void metGaussPiv(double A[n][n+1]){
	double z;
	int a;
	double x[n];
	//Llenado de la matriz
	int i;
	for(i=0;i<n;i++){
		for(int j=0;j<n+1;j++){
			cout<<i<<";"<<j<<": ";
			cin>>a;
			A[i][j]=a;
		}
		cout<<endl;
	}
	//imprimir
	imprimir(A);
	//Escalonado
	for(int j=0;j<=n;j++){
		for(i=0;i<=n-1;i++){
			if(i>j){
				z=A[i][j]/A[j][j];
				cout<<"Z: "<<z<<endl;
				for(int k=0;k<=n;k++)
					A[i][k]=int(A[i][k]-z*A[j][k]);

			}
		}
	}
	//imprimir
	imprimir(A);
	//sustitucion regresiva
	for(i=n-1;i>=0;i--){
		a=0;
		for(int k=i+1;k<=n-1;k++){
			a=a+A[i][k]*x[k];
		}
		x[i]=(A[i][n]-a)/A[i][i];
	}
	//imprimir var
	for(i=0;i<=n-1;i++)
        cout<<i+1<<" : "<<x[i]<<endl;
}
int main(){
    double A[n][n+1];
    metGaussPiv(A);

    return 0;
}
