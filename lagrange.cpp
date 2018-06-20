#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include "math.h"
#include <stdlib.h>
using namespace std;
ofstream archivo;
#define pi 3.141592654
#define ep 2.718281828
float x[50],y[50],z,l, valor=0;

void llenar_poli(float x[],float y[],int n){
	cout<<"ingrese x: "<<endl;
	for(int i=0; i<n; i++){
 		cin>>x[i];
    }
 
      cout<<endl;
      cout<<"ingrese y: "<<endl;
	for(int i=0; i<n; i++){
        cin>>y[i]; 
	}
}
void inter_lagrange(){
	int n,m;
	cout<<"ingrese el numero de elementos: "; cin>>n;
	llenar_poli(x,y,n);
	cout<<endl;
	cout<<"ingrese el numero de puntos a evaluar";
	cin>>m;cout<<endl;
	for(int i=0;i<m;i++){
		cin>>z;
		for(int i=0; i<n ;i++){
		    l=y[i;]
		    for(int j=0; j<n; j++){
		        if(i!=j){
		            l=(l*(z-x[j]))/(x[i]-x[j]);
		        }
		    }
		    valor=valor+l;
		}
		cout<<endl<<"Z= "<<z <<" es : "<<valor<<endl;
	}
}

int main(){
	inter_lagrange();
	return 0;
 }
