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
double A[n][n+1];
double L[n][n+1];
double U[n][n+1];
double x[n];
double y[n];
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
void imprimir(double A[n][n]){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
		cout<<endl;
	}
}

void llenar_matriz(double A[n][n+1]){
    int i;
    int a;
	for(i=0;i<n;i++){
		for(int j=0;j<n+1;j++){
			cout<<i<<";"<<j<<": ";
			cin>>a;
			A[i][j]=a;
		}
		cout<<endl;
	}
}
void escalona(double A[n][n+1]){
    int i;double z;
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
}

void identidad(double A[n][n+1]){
    int i;
	for(i=0;i<n;i++){
		for(int j=0;j<n+1;j++){
            A[i][j]=0;
		}
		cout<<endl;
	}
    for(i=0;i<n;i++){
		for(int j=0;j<n+1;j++){
            A[i][i]=1;
		}
		cout<<endl;
	}
}
void copia(double U[n][n+1],double A[n][n+1]){
    for(int i=0;i<n;i++){
        for(int j=0;j<n+1;j++){
            U[i][j]=A[i][j];
        }
    }
}
void susti_progresiva(double L[n][n+1]){
    int s=0;
    for(int i=0;i<n;i++){
        s=0;
        for(int k=0;k<i-1;k++){
            s=s+A[i][k]*y[k];
        }
        y[i]=(A[i][n]-s)/L[i][i];
    }
}

void llenarLU(){
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++)
            if(i>j){
              U[i][j]=0;    //Ceros debajo de la diagonal para la matriz U
            }
            else if(i==j){
              L[i][j]=1;    //Unos en la diagonal de L
            } else{
              L[i][j]=0;    //Ceros encima de la diagonal para la matriz L
            }
    }
}
void descomp_LU(double A[n][n+1]){
    int i,j,k,sum;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            sum=0;
            if(i<=j){               //Es decir, solo se llenan los elementos de la diagonal y encima de ella para U
                for(k=0;k<n;k++){
                    if(k!=i){
                        sum=sum+L[i][k]*U[k][j];
                    }else {
                        U[i][j]=A[i][j]-sum;               //Llenado de elementos restantes de U
                    }
                }
            } else{                 //Es decir, solo se llenan los elementos debajo de la diagonal para L
                for(k=0;k<n;k++){
                    if(k!=j){
                        sum=sum+L[i][k]*U[k][j];
                    }else {
                        L[i][j]=(A[i][j]-sum)/U[j][j];      //Llenado de elementos restantes de L
                    }
                }
            }
        }
    }

}
void susti_regresiva(double A[n][n+1]){
    int i,a;
    	double x[n];
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

void metDescompLU(double A[n][n+1]){
    llenarLU();
    llenar_matriz(A);
    descomp_LU(A);
    imprimir(L);
    imprimir(U);
    susti_progresiva(L);
    susti_regresiva(U);

}

void metGaussPiv(double A[n][n+1]){
    //Llenado de la matriz
    llenar_matriz(A);
	//imprimir
	imprimir(A);
	//Escalonado
	escalona(A);
	//imprimir
	imprimir(A);
	//sustitucion regresiva
	susti_regresiva(A);

}
int main(){
    //double A[n][n+1];
    //metGaussPiv(A);
    metDescompLU(A);


    return 0;
}
