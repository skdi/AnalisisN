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

/*int** escalonado(int **A,int n){

}*/
/*void imprimir(int &A[][],int n){
	for(i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
		cout<<endl;
	}
}*/
int main(){
	double z;int n=2;
	double A[n][n+1];
	int a;
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
	for(i=0;i<n;i++){
		for(int j=0;j<n+1;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
		cout<<endl;
	}
	//Escalonado
	for(int j=0;j<=n;j++){
		for(i=0;i<=n-1;i++){
			if(i>j){
				z=A[i][j]/A[j][j];
				cout<<"Z: "<<z<<endl;
				for(int k=0;k<=n;k++)
					A[i][k]=A[i][k]-z*A[j][k];
				
			}
		}
	}
	//imprimir
	for(i=0;i<n;i++){
		for(int j=0;j<n+1;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
		cout<<endl;
	}
	//sustitucion regresiva
	/*for(i=n;i>-1;i--){
		a=0;
		for(k=i+1;k<n){
			a=a+A[i][k]*x(k);
		}
		x[]
	}*/
	

}
