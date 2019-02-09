#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int a[100000+5];
ll sum[100000+5];
int main()
{
	int n,k,m;
	scanf("%d%d%d",&n,&k,&m);
	for (int i=1;i<=n;i++)
		scanf("%d",&a[i]);
	sort(a+1,a+n+1);
	sum[0]=0;
	for (int i=1;i<=n;i++)
		sum[i]=sum[i-1]+a[i];
	double ans=0;
	for (int i=0;i<min(m+1,n);i++)
	{
		ll tmp=sum[n]-sum[i];
		tmp+=min(1ll*(n-i)*k,1ll*(m-i));
		ans=max(ans,(double)tmp/(double)(n-i));
	}
	printf("%.10lf\n",ans);
	return 0;
}
