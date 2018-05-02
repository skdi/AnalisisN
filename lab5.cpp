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
double A[n][n];
double L[n][n];
double U[n][n];
double b[n];
double x[n];
double y[n];
void imprimir(double A[n][n]){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
		cout<<"B: "<<i<<": ";
		cout<<b[i]<<endl;
		//cout<<endl;
	}
}

void llenar_matriz(double A[n][n],double b[n]){
    int i;
    int a;
	for(i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<i<<";"<<j<<": ";
			cin>>a;
			A[i][j]=a;
		}
		cout<<i<<": ";
		cin>>a;
		b[i]=a;
		cout<<endl;
	}

}
void llenarLU(double L[n][n],double U[n][n]){
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
void descomp_LU(double A[n][n]){
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
void susti_regresiva(double A[n][n],double y[n]){
    int i,a;
    for(i=n-1;i>=0;i--){
		a=0;
		for(int k=i+1;k<=n-1;k++){
			a=a+A[i][k]*x[k];
		}
		x[i]=(y[i]-a)/A[i][i];
	}
	cout<<"Respuesta en X"<<endl;
	for(i=0;i<=n-1;i++)
        cout<<i+1<<" : "<<x[i]<<endl;
}
void susti_progresiva(double L[n][n],double b[n]){
    int s=0;
    for(int i=0;i<n;i++){
        s=0;
        for(int k=0;k<i-1;k++){
            s=s+L[i][k]*y[k];
        }
        y[i]=(b[i]-s)/L[i][i];
    }
    cout<<"Respuesta en Y"<<endl;
	for(int i=0;i<=n-1;i++)
        cout<<i+1<<" : "<<y[i]<<endl;
}

bool desc(double U[n][n]){
	bool a=1;
	for(int i=0;i<n;i++){

            	if(A[i][i]==0){
        			a=0;
			}  

    }
    return a;

}
void multiplicacion(double A[n][n],double B[n][n]){
	double C[n][n];
	for (int i=0;i<n;i++){
		for (int j=0;j<n;j++){
			 C[i][j]=0;
	         for(int k=0;k<n;k++){
	          C[i][j]=C[i][j]+A[i][k]*B[k][j];
	          }
	          cout<<"("<<i<<";"<<j<<"): "<<C[i][j]<<" ";
	       }
	       cout<<endl;

    }

}
void ResuelveSistConLU(double A[n][n],double b[n],double L[n][n],double U[n][n]){
    llenarLU(L,U);
    llenar_matriz(A,b);
    cout<<"MATRIZ A"<<endl;
	imprimir(A);
	descomp_LU(A);
    if(desc(U)){
	    cout<<"MATRIZ L"<<endl;
	    imprimir(L);
	    cout<<"MATRIZ U"<<endl;
	    imprimir(U);
	    susti_progresiva(L,b);
	    susti_regresiva(U,y);
	    cout<<"multiplicacion LU"<<endl;
	    multiplicacion(L,U);
    }
    else
    	cout<<"no tiene descomposicion LU"<<endl;





}
int main(){
	ResuelveSistConLU(A,b,L,U);

	return 0;
}
