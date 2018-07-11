#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
//funcion para el calculo de las constantes a0
float Aitmen(int i0, int i1, float* xAt, float* yAt) {
	if(i1-i0 == 0) {//si el punto a comparar es el mismo con el ya obtenido
    	return yAt[i0];//retornamos el punto
	}
	else {//sino
    	return (Aitmen(i0+1, i1, xAt, yAt) - Aitmen(i0, i1-1, xAt, yAt))/(xAt[i1]-xAt[i0]);//llamamos recursivamente a sus predecesores
	}
}
//funcion para interpolar los puntos x,y
float inter_new(float xAprox, int n, float* x, float* y) {
	float yAprox;
	float aitmenn[n];
	//calculamos las constantes a0...an
	for(int i = 0; i < n; i++) {
    	aitmenn[i] = Aitmen(0, i, x, y);
	}
	yAprox = 0;    
	for(int i = 0; i < n; i++) {
    	float aprox = aitmenn[i];//trabajamos la constante ai
    	for(int j = 0; j < i; j++) {//aproximacion
        	aprox *= (xAprox - x[j]);
    	}
   	 
    	yAprox += aprox;
	}
    
	return yAprox;
}
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
void main_newton(){
	float xAprox, yAprox;
	int n,m;
	cout << "ingrese el numero de elementos: ";
	cin >> n;
	float x[n];
	float y[n];
	llenar_poli(x,y,n);
	cout << "ingrese el numero de puntos a evaluar ";
	cin>>m;
	for(int i=0;i<m;i++){
    	cin >> xAprox;
    	yAprox = inter_new(xAprox, n, x, y);
    	cout << "f(x) = " << yAprox<<endl;
	}

}
int main() {
	main_newton();

    
}


