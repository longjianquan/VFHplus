#include <iostream>
#include <eigen3/Eigen/Dense>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <math.h>
#include <fstream>
#define PI 3.1415926
const int step=1;
const int f=5;
const int dmax=10;
const int smax=18;
const float b=2.5;
const float a=1+b*(pow(dmax,2));
const int C=15;
const double alpha=f*PI/180.0;
double n=360/f;
double threshold=1000.0;
const double rsafe=0.3;
double kb=90/f;
using namespace std;
using namespace Eigen;
struct pointpos
{
    Vector2f pos;
};

template<class T>
class VFHplus
{
public:
    VFHplus();
    vector<vector<T> >Readmap();
    vector<pointpos>GetObstacle(vector<vector<T> >&map);
    vector<pointpos> Solution(vector<vector<T> >&map,vector<pointpos>&obstacle);
    double caculatebeta(Vector2f &robot,Vector2f &goal);
    double caculatedistance(Vector2f &a,Vector2f &b);
    double rad2deg(double &f);
    double howmany(double &c1,double &c2);
};
template<class T>
VFHplus<T>::VFHplus()
{
    cout<<"start"<<endl;
}
template<class T>
vector<vector<T> > VFHplus<T>::Readmap()
{
    ifstream infile;
    infile.open("/home/ljq/test/1/1.txt");
    while (infile.good()&& !infile.eof())
    {
       vector< vector<int> >map;
       pointpos pp;
       for(int i=0;i<15;i++)
       {
           vector<int>temp;
           temp.clear();
           for(int j=0;j<15;j++)
           {
               int tempvalue;
               infile>>tempvalue;
               temp.push_back(tempvalue);
           }
           map.push_back(temp);
       }
        return map;
    }
  infile.close();
}
template<class T>
vector<pointpos> VFHplus<T>::GetObstacle(vector<vector<T> >&map)
{
    vector<pointpos>obstacle;
    pointpos pp;
    for(int i=0;i<map.size();i++)
    {
        for(int j=0;j<map[0].size();j++)
        {
            if(map[i][j]==1)
            {
                pp.pos.x()=i;
                pp.pos.y()=j;
                obstacle.push_back(pp);
            }
        }
    }
    return obstacle;
}
template<class T>
double VFHplus<T>::caculatebeta(Vector2f &robot,Vector2f &goal)
{
    double dy=goal.y()-robot.y();
    double dx=goal.x()-robot.x();
    double beta;
    if (dx==0)
    {
        beta=PI/2;
    }
    else
    {
        beta=atan(dy/dx);
        if(dx<0)
        {
            if(dy>0)
            {
                beta=PI-abs(beta);
            }
            else
            {
                beta=PI+abs(beta);
            }
        }
        else
        {
            if(dy<0)
            {
                beta=2*PI-abs(beta);
            }
        }
    }

    return beta;
}
template<class T>
double VFHplus<T>::caculatedistance(Vector2f &a,Vector2f &b)
{
    double dist=sqrt(pow(a.x()-b.x(),2)+pow(a.y()-b.y(),2));
    return dist;
}
template<class T>
double VFHplus<T>::rad2deg(double &f)
{
    double value=f*180/PI;
    return value;
}
template<class T>
double VFHplus<T>::howmany(double &c1, double &c2)
{
    double dif=min(min(abs(c1-c2),abs(c1-c2-n)),abs(c1-c2+n));
    return dif;
}
template<class T>
vector<pointpos> VFHplus<T>::Solution(vector<vector<T> >&map,vector<pointpos>&obstacle)
{
    Vector2f start(0.0,0.0);
    Vector2f goal(14.0,14.0);
    Vector2f robot;
    robot=start;
    double kt=round(caculatebeta(robot,goal)/alpha);
   //cout<<kt<<endl;

    if(kt==0.00)
    {
        kt=n;
    }
    while (caculatedistance(robot,goal)!=0.00)
    {
        if (caculatedistance(robot,goal)>step*3)
        {
            int i=0;
            vector<double>mag(n,0);
            vector<double>his(n,0);
            while (i<obstacle.size())
            {
                Vector2f temp;
                temp.x()=obstacle[i].pos.x();
                temp.y()=obstacle[i].pos.y();
                double d=caculatedistance(robot,temp);
                if(d<dmax)
                {
                    double beta=caculatebeta(robot,temp);
                    double rangle=asin(rsafe/d);
                    double k=round(beta/alpha);
                    if (k==0.00)
                    {
                        k=n;
                    }
                    vector<double>h(n+2,0);
                    if((5*k>rad2deg(beta)-rad2deg(rangle))&&(5*k<rad2deg(beta)+rad2deg(rangle)))
                    {
                        h[k]=1;
                    }
                    else
                    {
                        h[k]=0;
                    }
                    double m=pow(C,2)*(a-b*(pow(d,2)));
                    mag[k]=max(mag[k],m*h[k]);
                    i=i+1;
                }
                else
                    i=i+1;
            }
            his=mag;
            int j=0,q=0;
            vector<double>c(n+2,0);
            while(q<n)
            {
                if(his[q]<threshold)
                {
                    double kr=q;
                    double kl;
                    while(q<=n && his[q]<threshold)
                    {
                        kl=q;
                        q=q+1;
                    }

                    if(kl-kr>smax)
                    {
                        c[j]=round(kl-smax/2);
                        c[j+1]=round(kr+smax/2);
                        j=j+2;
                        if((kt>=kr)&&(kt<=kl))
                        {
                            c[j]=kt;
                            j=j+1;
                        }
                    }
                    else if (kl-kr>smax/5)
                    {
                        c[j]=round((kr+kl)/2);
                        j=j+1;
                    }
                }
                else
                {
                    q=q+1;
                }
            }
            vector<double>g(j,0);
            vector<double>how(j,0);
            int minfind=9999999;
            int ft=0;
            for(int i=0;i<j-1;i++)
            {
                g[i]=c[i];
                how[i]=5*howmany(g[i],kt)+2*howmany(g[i],kb)+2*howmany(g[i],kb);
                if(how[i]<minfind)
                {
                    minfind=how[i];
                    ft=i;
                }
            }
            int fk=ft;
            double kb=g[fk];
            double ccc=step*cos(kb*alpha);
            double ddd=step*sin(kb*alpha);

            robot.x()=robot.x()+step*cos(kb*alpha);
            robot.y()=robot.y()+step*sin(kb*alpha);
            cout<<robot.x()<<" "<<robot.y()<<endl;
            double kt=round(caculatebeta(robot,goal)/alpha);
            if(kt==0)
            {
                double kt=n;
            }

        }
        else
        {
            break;
        }
    }

}
int main()
{ 
    VFHplus<int>VFHplus;

    vector< vector<int> >map=VFHplus.Readmap();
    vector<pointpos>obstacles=VFHplus.GetObstacle(map);
    VFHplus.Solution(map,obstacles);

}
