#include <bits/stdc++.h>
using namespace std;
int main()
{
	int n,r;
	scanf("%d%d",&n,&r);
	double k=2*acos(-1)/n;
	double a=(-1-cos(k)),b=r*(2-2*cos(k)),c=r*r*(1-cos(k));
	double ans=(-b-sqrt(b*b-4*a*c))/(2*a);
	printf("%.10lf\n",ans);
	return 0;
}
