#include <iostream>
#include <math.h>
#include <cmath>
//#define _USE_MATH_DEFINES
#define pipi 3.14159265358979
using namespace std;
double f1(double x)
{
    return pow(x,4);
    //integral x^4 dx de [0,1] = 1/5
}
double f2(double x)
{
    return pow(x,5);
}
double f(double x)
{
    return sin(x);
    //integral sen(x)dx de [0,1] = 1
}
double Trapecio(double a,double b,double( *f)(double),int n)
{

    double h=(b-a)/n;
    double suma=0,r;
    for(int i=1;i<=n-1;i++)
    {
        suma=suma+(2*f(a+(i*h)));
    }
    r=(h/2)*(f(a)+suma+f(b));
    return r;
}
double Simpson13(double a,double b,double(*f)(double),int n)
{
    if(n%2==0)
    {
        double h=(b-a)/n;
        int xjpar=0,xiimpar=0;
        double s1=0,s2=0,rs=0,rs1;
        int limite1=(n-2)/2;

        if(limite1!=0)
        {
            for(int i=1;i<=n/2;i++)
            {
                xiimpar=2*i-1;
                //cout<<xiimpar<<" ";
                s2=s2+f(a+h*xiimpar);
            }
            //cout<<"S2:"<<s2;
            //cout<<endl;
            for(int j=1;j<=((n/2)-1);j++)
            {
                xjpar=2*j;
                //cout<<xjpar<<" ";
                s1=s1+f(a+h*xjpar);
            }
            //cout<<"S1:"<<s1;

            rs=(h/3)*(f(a)+2*s1+4*s2+f(b));
            return rs;
        }
        else
        {
            int d;
            s1=0;
            for(int k=1;k<=n/2;k++)
            {
                d=2*k-1;
                s2=s2+f(a+h*d);
            }
            rs1=(h/3)*(f(a)+f(b)+4*s2);
            return rs1;
        }
    }
    else
    {
        cout<<" 'n' no es par "<<endl;
    }
}
double Simpson38(double a,double b,double (*f)(double),int n)
{
    if(n%3==0)
    {
        double h=(b-a)/n;
        double xi=0;
        double s1=0,s2=0,s3=0,rs;
        int limite2=n-3;

        if(limite2!=0)
        {
            int i=1,j=2,k=3;
            while(k<=n-3)
            {
                s3=s3+f(a+h*k);
                //cout<<k<<" ";
                k=k+3;
                while(i<=n-2)
                {
                    s1=s1+f(a+h*i);
                    //cout<<i<<" ";
                    i=i+3;
                    while(j<=n-1)
                    {
                        s2=s2+f(a+h*j);
                        //cout<<j<<" ";
                        j=j+3;
                    }
                }
            }
        }
        else
        {
            s3=0;
            int x=1,y=2;
            while(y<=n-1)
            {
                s2=s2+f(a+h*y);
                y=y+3;
                while(x<=n-2)
                {
                    s1=s1+f(a+h*x);
                    x=x+3;
                }
            }
        }
        rs=((3*h)/8)*(f(a)+3*s1+3*s2+2*s3+f(b));
        return rs;
    }
    else
    {
        //cout<<"'n' no es multiplo de 2";
    }
}
void viewTable(int n,double &a,double &b)
{
    cout<<endl;
    cout<<"Funcion: sen(x)"<<endl;
    cout<<"Intervalo: "<<"["<<a<<";"<<b<<"]";
    cout<<endl;
    cout<<"\n";
    cout<<" n"<<'\t'<<"Trapecio"<<'\t'<<"Simpson 1/3"<<'\t'<<"Simpson 3/8"<<endl;
    cout<<"\n";
    cout<<"-------------------------------------------------------------"<<endl;
    for(int i=6;i<193;i=i*2)
    {
        cout<<" "<<i<<"\t"<<Trapecio(a,b,f,i)<<"     \t"<<Simpson13(a,b,f,i)<<"    \t"<<Simpson38(a,b,f,i)<<endl;
        //cout<<i<<" ";
    }
}
int main()
{
    double a=0;
    double b=pipi/2;
    int n=6;
    cout<<"\n  Tabla de Comparacion entre las Reglas de: \n"<<endl;
    cout<<".........................................................."<<endl;
    viewTable(n,a,b);
    //cout<<Simpson13(a,b,f,n);
    return 0;
}
