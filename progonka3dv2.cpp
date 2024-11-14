#include <iostream>
#include <cmath>
using namespace std;

int main()
{
	double a[5]={0,-3,-5,-6,-5};
	double b[5]={2,8,12,18,10};
	double c[5]={-1,-1,2,-4,0};
	double d[5]={-25,72,-69,-156,20};
	int n=5;
	double j[n];
	double alfa[n];
	double betta[n];
	double x[n];
	cout<<"Straight"<<'\n';
	for (int i=0; i<n; i++)//straight
	{
		cout<<"Input:"<<'\n';
		cout<<"a["<<i<<"]="<< a[i]<<'\n';
		cout<<"b["<<i<<"]="<< b[i]<<'\n';
		cout<<"c["<<i<<"]="<< c[i]<<'\n';
		j[i]=b[i]+a[i]*alfa[i-1];
		alfa[i]=-c[i]/j[i];
		betta[i]=(d[i]-a[i]*betta[i-1])/j[i];
		cout<<"Output:"<<'\n';
		cout<<"j["<<i<<"]="<< j[i]<<'\n';
		cout<<"alfa["<<i<<"]="<< alfa[i]<<'\n';
		cout<<"betta["<<i<<"]="<< betta[i]<<'\n';
	}
	cout<<"----------------------------------------"<<'\n';
	cout<<"Back"<<'\n';
	for(int i=4; i>=0; i--)//back
	{
		cout<<"Input:"<<'\n';
		cout<<"x["<<i+1<<"]="<< x[i+1]<<'\n';
		cout<<"alfa["<<i<<"]="<< alfa[i]<<'\n';
		cout<<"betta["<<i<<"]="<< betta[i]<<'\n';
		cout<<"Output:"<<'\n';		
		x[i]=alfa[i]*x[i+1]+betta[i];
		cout<<"x["<<i<<"]="<< x[i]<<'\n';
	}
	
}
