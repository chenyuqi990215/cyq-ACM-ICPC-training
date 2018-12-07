#include <bits/stdc++.h>
using namespace std;
const int maxn=1e5+5;
int a[maxn];
int main()
{
	int n,m;
	scanf("%d%d",&n,&m);
	long long sum=0;
	for (int i=0;i<n;i++)
	{
		scanf("%d",&a[i]);
		sum+=a[i];
	}
	sort(a,a+n);
	reverse(a,a+n);
	long long ans=0;
	int h=a[0];
	for (int i=0;i<n;i++)
	{
		if (i==n-1) ans=ans+h;
		else
		{
			ans=ans+max(1,h-a[i+1]);
			h=max(1,h-max(1,h-a[i+1]));
		}
	}
	printf("%I64d\n",sum-ans);
	return 0;
}
