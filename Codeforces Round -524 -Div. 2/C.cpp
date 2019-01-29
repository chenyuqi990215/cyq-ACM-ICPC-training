#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll calb(ll x1,ll y1,ll x2,ll y2)
{
	ll w=x2-x1+1,h=y2-y1+1;
	if (w<=0 || h<=0) return 0;
	if ((w*h)%2==0) return w*h/2;
	if ((x1+y1)%2==0) return (w*h-1ll)/2;
	return (w*h+1ll)/2;
}
ll calw(ll x1,ll y1,ll x2,ll y2)
{
	ll w=x2-x1+1,h=y2-y1+1;
	if (w<=0 || h<=0) return 0;
	if ((w*h)%2==0) return w*h/2;
	if ((x1+y1)%2==0) return (w*h+1ll)/2;
	return (w*h-1ll)/2;
}
ll cals(ll x1,ll y1,ll x2,ll y2)
{
	ll w=x2-x1+1,h=y2-y1+1;
	if (w<=0 || h<=0) return 0;
	return w*h;
}
int main()
{
	int t,n,m,x1,x2,x3,x4,y1,y2,y3,y4;
	scanf("%d",&t);
	while (t--)
	{
		scanf("%d%d",&n,&m);
		scanf("%d%d%d%d",&x1,&y1,&x2,&y2);
		scanf("%d%d%d%d",&x3,&y3,&x4,&y4);
		ll totw=0;
		totw+=calw(1,1,m,n);
		totw+=calb(x1,y1,x2,y2);
		totw-=calw(x3,y3,x4,y4);
		ll x5=max(x1,x3),y5=max(y1,y3),x6=min(x2,x4),y6=min(y2,y4);
		totw-=cals(x5,y5,x6,y6);
		totw+=calw(x5,y5,x6,y6);
		ll totb=cals(1,1,m,n)-totw;
		printf("%lld %lld\n",totw,totb);
	}
	return 0;
}
