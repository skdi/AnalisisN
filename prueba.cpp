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


/*BISECCION LAB! y FALSA POSICION LAB2*/
double fun( double x){
    return ((3*(pow(x,5))+2*x*(pow(ep,-x))+5*cos(x))-4);
    }

double fun2(double x){
    return (4*x-5+pow(ep,-2*x+1)*sin(x+pi)-pow(x,2));
    }

double fun3(double x){
    return (10*pow(x,3)+x*atan(x+2*pi)+5);
}

double fun4(double x){
    return (pow(x,2)*cos(x*2-1)+pow(ep,x-2)-x);
}


double biseccion_iterador(double funcion(double X), double parametro1, double parametro2, double iteracion){
    archivo.open("datos.txt");
    archivo<< "i" << "\t" << "parametro1" << "\t" << "parametro2" << " \t" << "iteracion"<< " \t" << "funcion"<<endl;
    double m = (parametro1+parametro2)/2;
    for(int i =1; i<=iteracion; i++){
        archivo<< i << "\t  " << parametro1 << "\t  " << parametro2 << "\t  " << m << "\t   " << funcion(m)<<endl;

        if (funcion(parametro1)*funcion(m)<0){
            parametro2=m;
        }
        else{
            parametro1=m;
        }
        m=(parametro1+parametro2)/2;

    }
    return m;

}

double biseccion_tolerancia(double funcion(double X), double parametro1, double parametro2, double tolerancia){
    archivo.open("datos.txt");
    archivo<< "i" << "\t" << "parametro1" << "\t" << "parametro2" << " \t" << "iteracion"<< " \t" << "funcion"<<endl;
    double m = (parametro1+parametro2)/2;
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
        m = (parametro1+parametro2)/2;

    }
    return m;
}

double biseccion_tolerancia2(double funcion(double X), double parametro1, double parametro2, double tolerancia){
    archivo.open("datos.txt");
    archivo<< "i" << "\t" << "parametro1" << "\t" << "parametro2" << " \t" << "iteracion"<< " \t" << "funcion"<<endl;
    double m = (parametro1+parametro2)/2;
    int i=0;
    while(abs(parametro1-parametro2)>tolerancia){
        archivo<< i << " " << parametro1 << " " << parametro2 << " " << tolerancia << " " << funcion(tolerancia)<<endl;
        if (funcion(parametro1)*funcion(m)<0){
            parametro2=m;
        }
        else{
            parametro1=m;
        }
        i++;
        m = (parametro1+parametro2)/2;

    }
    return m;
}

double falsaposicion_iterador(double funcion(double X), double parametro1, double parametro2, double iteracion){
    archivo.open("datos.txt");
    archivo<< "i" << "\t" << "parametro1" << "\t" << "parametro2" << " \t" << "iteracion"<< " \t" << "funcion"<<endl;
    double m = parametro2-funcion(parametro2)*((parametro2-parametro1)/(funcion(parametro2)-funcion(parametro1)));
    for(int i =1; i<=iteracion; i++){
        archivo<< i << "\t  " << parametro1 << "\t  " << parametro2 << "\t  " << m << "\t   " << funcion(m)<<endl;

        if (funcion(parametro1)*funcion(m)<0){
            parametro2=m;
        }
        else{
            parametro1=m;
        }
        m = parametro2-funcion(parametro2)*((parametro2-parametro1)/(funcion(parametro2)-funcion(parametro1)));

    }
    return m;

}
double falsaposicion_tolerancia(double funcion(double X), double parametro1, double parametro2, double tolerancia){
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
}

double falsaposicion_tolerancia2(double funcion(double X), double parametro1, double parametro2, double tolerancia){
    archivo.open("datos.txt");
    archivo<< "i" << "\t" << "parametro1" << "\t" << "parametro2" << " \t" << "iteracion"<< " \t" << "funcion"<<endl;
    double m = parametro2-funcion(parametro2)*((parametro2-parametro1)/(funcion(parametro2)-funcion(parametro1)));
    int i=0;
    while((parametro2-parametro1)>tolerancia){
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
}


void main_biseccion_falsaposicion(){

    //BISECCION
    //FUNCION 1
    //Metodo 1
    //biseccion_iterador(fun,-2,2,100);

    //Metodo 2
    //biseccion_tolerancia(fun,-2,2,0.0000001);

    //Metodo 3
    //biseccion_tolerancia2(fun,-2,2,0.0000001);


        //FUNCION 2
    //Metodo 1
    //biseccion_iterador(fun2,-4,-3,100);

    //Metodo 2
    //biseccion_tolerancia(fun2,-4,-3,0.0000001);

    //Metodo 3
    //biseccion_tolerancia2(fun2,-4,-3,0.0000001);
////////////////////////////////////////////////////////////////////////////////////
    //    FALSA POSICION
    //FUNCION 1
    //Metodo 1
    //falsaposicion_iterador(fun3,-1,1,100);

    //Metodo 2
    //falsaposicion_tolerancia(fun3,-1,1,0.0000001);

    //Metodo 3
    //falsaposicion_tolerancia2(fun3,-1,1,0.0000001);


        //FUNCION 2
    //Metodo 1
    //falsaposicion_iterador(fun4,0,1,100);

    //Metodo 2
    //falsaposicion_tolerancia(fun4,0,1,0.0000001);

    //Metodo 3
    falsaposicion_tolerancia(fun4,0,1,0.0001);


}

/*Newton lab3*/

double funN1( double x){
    return (1/(pow(3,x)));
    }
double funN2(double x){
    return (2/pow(3,pow(2,x)));
}
double funN3(double x){
    return (1/pow(x,2));
}
double funN4(double x){
    if(x>=0){
        if(x==0)
            return 2;
        else{
            return ((1/2*funN4(x-1))+(1/funN4(x-1)));
        }
}
}
double funN5(double x){
    return pow(x,3)+4*pow(x,2)-10;
}
double d_funN5(double x){
    float h=0.0001;
    return ((funN5(x+h)-funN1(x))/h);
}
double newton(double r,double a,double b){
    if(r>=a && r<=b)
        return r-(funN5(r)/d_funN5(r));
}
double alpha_1(int n){
    return ((log(abs(funN1(n-1)))-log(abs(funN1(n))))/(log(abs(funN1(n)))/log(abs(funN1(n-1)))));
}
double alpha_2(int n,double x0=0){
    return ((log(abs(funN2(n-1)-x0))-log(abs(funN2(n)-x0)))/(log(abs(funN2(n)-x0))/log(abs(funN2(n-1)-x0))));
}
double alpha_3(int n,double x0=0){
    return ((log(abs(funN3(n-1)-x0))-log(abs(funN3(n)-x0)))/(log(abs(funN3(n)-x0))/log(abs(funN3(n-1)-x0))));
}
double alpha_4(int n,double x0=0){
    return ((log(abs(funN4(n-1)-x0))-log(abs(funN4(n)-x0)))/(log(abs(funN4(n)-x0))/log(abs(funN4(n-1)-x0))));
}

void main_newton(){
    double i=1;
    int n=2;
    double r0;
    double r=1.5;
    while(n<100){
        i=abs(alpha_4(n)-alpha_4(n-1));
        if(i>0.001){
            cout<<n-1<<":"<<alpha_4(n)<<" Error: "<<i<<endl;
            n++;
        }   else break;
        //NEWTON
        r0=r;
        n++;
        r=newton(r,1,2);
        i=abs(r-r0);
        cout<<n-1<<": "<<r<<" Error: "<<i<<endl;



    }

}

/*Gauss-lab4*/

#define n 2

void imprimir(double A[n][n+1]){
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
        }
        cout<<endl;
    }
}

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

    //permutacion
    int fmay;double may;double paso;
    for(int i=0;i<=n;i++){
        may=abs(A[i][i]);
        fmay=i;
        for(int j=i;j<=n;j++){
            if(abs(A[j][i])>may){
                may=abs(A[j][i]);
                fmay=j;
            }
        }
        for(int k=0;k<=n;k++){
            paso=A[i][k];
            A[i][k]=A[fmay][k];
            A[fmay][k]=paso;
        }

    }
    //Escalonado
    double temp;
    for(int j=0;j<=n;j++){
        for(i=0;i<=n-1;i++){
            if(i>j){
                z=A[i][j]/A[j][j];
                cout<<"Z: "<<z<<endl;
                for(int k=0;k<=n;k++)
                    if(A[i][i]==0)
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
void main_gauss_lab4(){
    double A[n][n+1];
    metGaussPiv(A);


}

/*Descomposicion LU*/

double A[n][n];
double L[n][n];
double U[n][n];
double b[n];
double x[n];
double y[n];

void permutar(double A[n][n],double b[n]){
    int fmay;double may;double paso;double temp;
    for(int i=0;i<n;i++){
        may=abs(A[i][i]);
        fmay=i;
        for(int j=i;j<n;j++){
            if(abs(A[j][i])>may){
                may=abs(A[j][i]);
                fmay=j;
            }
        }
        for(int k=0;k<n;k++){
            paso=A[i][k];
            A[i][k]=A[fmay][k];
            A[fmay][k]=paso;
            temp=b[i];
            b[i]=b[fmay];
            b[fmay]=temp;
        }


    
}}
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
    cout<<"Matriz A permutada"<<endl;
    permutar(A,b);
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
int main_LU(){
    ResuelveSistConLU(A,b,L,U);

    return 0;
}

/*Jacobi*/

int dim;
float norma(float vector1[],float vector2[]);
float suma_jacobi(float Matriz[], float vector[], int componente);

int main_jacobi(){
    int i,j,iteraciones=0;
    float error,epsilon;
    cout<<"METODO DE JACOBI DE RESOLUCION DE SISTEMAS Ax=b "<<endl;

    cout<<"Dimension de la matriz A: ";
    cin>>dim;
    float A[dim][dim],b[dim],x[dim],x_prev[dim],aux[dim];

    cout<<endl<<" Elementos de la matriz A: "<<endl;
    for(i=0;i<dim;i++) for(j=0;j<dim;j++){
        cout<<"A("<<i<<","<<j<<")"; cin>>A[i][j]; cout<<endl;
    }

    cout<<" Elementos del vector b: "<<endl;
    for(i=0;i<dim;i++){
        cout<<"b("<<i<<")="; cin>>b[i];
    }

    cout<<"\n Error de parada: \n";
    cout<<"E= "; cin>>epsilon;
    error=epsilon+1;

    //cominezo algoritmo de Jacobi
    //Error se mide como la norma del vector diferenceia entre la iteracion i e i+1
    cout<<"\n Valor inicial de la iteracion: "<<endl;
    for(i=0;i<dim;i++){
        cout<<"x0("<<i<<")="; cin>>x_prev[i];
    }
    while (error>epsilon){
        for(i=0;i<dim;i++){
            for(j=0;j<dim;j++) 
                aux[j]=A[i][j];//copia la fila a trabajar
                x[i]=(1/A[i][i])*(b[i]-suma_jacobi(aux,x_prev,i));
        }
        error=norma(x,x_prev);

        cout<<"\n\n Iteracion"<<iteraciones<<endl;
        for(i=0;i<dim;i++){
            x_prev[i]=x[i];
            cout<<"X("<<i<<")="<<x[i]<<endl;
        }

        iteraciones++;
        if (iteraciones==10) error=epsilon-1;
    }

    cout<<"Solucion del sistema"<<endl;
    cout<<"Numero de iteraciones: "<<iteraciones<<endl;
    for(i=0;i<dim;i++){
        cout<<"x("<<i<<")="<<x[i]<<endl;
    }
    return 1;
}
float norma2(float vector1[],float vector2[]){
    float aux=0;
    for(int i=0;i<dim;i++){
        aux+=abs(vector1[i]-vector2[i]);
    }
    return aux/3;

}    
float norma(float vector1[],float vector2[]){
    float aux=0;
    int i;
    for(i=0;i<dim;i++){
        aux=aux+(vector1[i]-vector2[i])*(vector1[i]-vector2[i]);
    }
    return aux;
}

float suma_jacobi(float Matriz[], float vector[], int componente)
{
    float aux=0;
    int i;
    for(i=0;i<dim;i++){
        if (componente!=i){
            aux=aux+Matriz[i]*vector[i];
        }
    }
    return aux;
}

/*Gauss seidel*/
int main_gauss_seidel(){
    int dim,i,j,k,flag=0;int count,l,e;
    double y,temp,temp2;    float error,epsilon;
    cout<<"Dimension de la matriz A: ";
    cin>>dim;
    float A[dim][dim],b[dim],x[dim],x_prev[dim],aux[dim];

    cout<<endl<<" Elementos de la matriz A: "<<endl;
    for(i=0;i<dim;i++) for(j=0;j<dim;j++){
        cout<<"A("<<i<<","<<j<<")"; cin>>A[i][j]; cout<<endl;
    }

    cout<<" Elementos del vector b: "<<endl;
    for(i=0;i<dim;i++){
        cout<<"b("<<i<<")="; cin>>b[i];
    }

    cout<<"\n Error de parada: \n";
    cout<<"E= "; cin>>epsilon;
    error=epsilon+1;

    //cominezo algoritmo de Jacobi
    //Error se mide como la norma del vector diferenceia entre la iteracion i e i+1
    cout<<"\n Valor inicial de la iteracion: "<<endl;
    for(i=0;i<dim;i++){
        cout<<"x0("<<i<<")="; cin>>x[i];
    }
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
        }
        cout<<"B: "<<i<<": ";
        cout<<b[i];
        cout<<" X: "<<i<<": ";
        cout<<x[i]<<endl;
    }
    //diagonalizar la matriz
    for(i=0;i<dim;i++){
        for(k=i+1;k<dim;k++){
            if(abs(A[i][i])<abs(A[k][i])){
                for(j=0;j<dim;j++){
                    temp=A[i][j];
                    A[i][j]=A[k][j];
                    A[k][j]=temp;
                    temp2=b[i];
                    b[i]=b[k];
                    b[k]=temp2;
                }
            }
        }
    }
    cout<<"A diagonalizada: ";
    for(i=0;i<dim;i++){
        for(j=0;j<dim;j++){
            cout<<"("<<i<<";"<<j<<"): "<<A[i][j]<<" ";
        }
        cout<<"B: "<<i<<": ";
        cout<<b[i];
    }

    cout<<"iteracion N:";
    for(i=0;i<dim;i++)
        cout<<"x "<<i+1<<endl;
    do{
        cout<<count+1;
        for(i=0;i<dim;i++){
            y=x[i];
            x[i]=A[i][dim];
            for(j=0;j<dim;j++){
                if(j!=i)
                    x[i]=x[i]-A[i][j]*x[j];
            }
            x[i]=x[i]/A[i][i];
            if(abs(x[i]-y)<=epsilon)
                flag++;
            cout<<"|"<<x[i]<<"|";
        }
        cout<<endl;
        count++;

    }while(flag<dim);
    cout<<"La solucion es la sig"<<endl;
    for(i=0;i<dim;i++)
        cout<<"x"<<i+1<<": "<<x[i]<<endl;
}
int main(){
    cout<<"Metodos numericos: "<<endl;
    //main_biseccion_falsaposicion();
    //main_newton();
    //main_gauss_lab4();
    //main_LU();
    //main_jacobi();
    main_gauss_seidel();
    archivo.close();
    return 0;
    }
