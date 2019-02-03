#include <bits/stdc++.h>
#define mp make_pair
using namespace std;
typedef pair<int,int> pii;
typedef long long ll;
int a[100000+5],c[100000+5];
pii p[100000+5];
int n;
ll get(int &q,int d)
{
	ll ret=0;
	while (q<=n)
	{
		if (a[p[q].second]>=d)
		{
			ret+=1ll*d*c[p[q].second];
			a[p[q].second]-=d;
			return ret;
		}
		else
		{
			ret+=1ll*a[p[q].second]*c[p[q].second];
			d-=a[p[q].second];
			a[p[q].second]=0;
			q++;
		}
	}
	return 0;
}
int main()
{
	int m,t,d;
	scanf("%d%d",&n,&m);
	for (int i=1;i<=n;i++)
		scanf("%d",&a[i]);
	for (int i=1;i<=n;i++)
		scanf("%d",&c[i]);
	for (int i=1;i<=n;i++)
		p[i]=mp(c[i],i);
	sort(p+1,p+n+1);
	int q=1;
	for (int i=0;i<m;i++)
	{
		scanf("%d%d",&t,&d);
		ll ans=0;
		if (a[t]>=d)
		{
			ans+=1ll*c[t]*d;
			a[t]-=d;
		}
		else
		{
			ans+=1ll*c[t]*a[t];
			d-=a[t];
			a[t]=0;
			ll temp=get(q,d);
			if (temp==0)
				ans=0;
			else ans+=temp;
		}
		printf("%lld\n",ans);
	}
	return 0;
}
