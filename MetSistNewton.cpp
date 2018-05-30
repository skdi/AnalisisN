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
#define n 2
double A[n][n];
double L[n][n];
double U[n][n];
double J[n][n];
double F[n];
double b[n];
double x[n];
double y[n];
void imprimir(double A[n][n]){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
		}
        cout<<endl;
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
                         if(U[i][i]==0){
                            //permutacion
                            int fmay;double may;double paso;
                            for(int i=0;i<=n;i++){
                                may=abs(U[i][i]);
                                fmay=i;
                                for(int j=i;j<=n;j++){
                                    if(abs(U[j][i])>may){
                                        may=abs(U[j][i]);
                                        fmay=j;
                                    }
                                }
                                for(int k=0;k<=n;k++){
                                    paso=U[i][k];
                                    U[i][k]=U[fmay][k];
                                    U[fmay][k]=paso;
                                }
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
    /*cout<<"Respuesta en Y"<<endl;
	for(int i=0;i<=n-1;i++)
        cout<<i+1<<" : "<<y[i]<<endl;*/
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
	descomp_LU(A);
    if(1){
	    susti_progresiva(L,b);
	    susti_regresiva(U,y);
    }
    else
    	cout<<"no tiene descomposicion LU"<<endl;

}
//////////////////////////////////////////////////////////////////////////
//NEWTON
//////////////////////////////////////////////////////////////////////////

double f1(double x1,double x2){
    return pow(x1,2)+pow(x2,2)-2;
}
double f2(double x1,double x2){
    return x1-x2-1;
}
double df1x(double x,double y){
    float h=0.0001;
    return ((f1(x+h,y)-f1(x,y))/h);
}
double df1y(double x,double y){
    float h=0.0001;
    return ((f1(x,y+h)-f1(x,y))/h);
}
double df2x(double x,double y){
    float h=0.0001;
    return ((f2(x+h,y)-f2(x,y))/h);
}
double df2y(double x,double y){
    float h=0.0001;
    return ((f2(x,y+h)-f2(x,y))/h);
}
void F1(double x,double y){
    F[0]=f1(x,y);
    F[1]=f2(x,y);
}
void JF1(double x,double y){
    J[0][0]=df1x(x,y);
    J[0][1]=df1y(x,y);
    J[1][0]=df2x(x,y);
    J[1][1]=df2y(x,y);
}



void MetNewtonSist(double b[n],double x0[n],double tol){
    double error=1,temp1,temp2;
    double x1,x2;//x3
    x1=x0[0];x2=x0[1];//x3=x0[2];
    while(error>tol){
        temp1=temp2=0;
        x1=f1(x1,x2);x2=f2(x1,x2);//tenemos F(Xn-1)
        JF1(x1,x2);//tenemos J
        imprimir(J);
        ResuelveSistConLU(J,b,L,U);//tenemos  X(n)
        //NORMA
        for(int i=1;i<n;i++){
            temp1+=x[i];
            temp2+=x[i-1];
        }
        temp1/=n;temp2/=n;
        error=abs(temp1-temp2);

    }
    


}


int main(){
    double b[n],x0[n],tol;
    for(int i=0;i<n;i++){
        cout<<i<<": ";
        cin>>b[i];
        cout<<i<<": ";
        cin>>x0[i];
    }
    cout<<"ingrese la tolerancia: ";cin>>tol;
	MetNewtonSist(b,x0,tol);

	return 0;
}
