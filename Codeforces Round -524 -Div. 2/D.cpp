#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int main()
{
	int t,n;
	ll k;
	scanf("%d",&t);
	while (t--)
	{
		scanf("%d%lld",&n,&k);
		if (n>=32)
		{
			printf("YES %d\n",n-1);
		}
		else
		{
			bool flag=false;
			for (int a=0;a<=n;a++)
			{
				ll l=(1ll<<(n-a+1))-2-(n-a);
				ll r=((1ll<<(2*n))-1)/3-((1ll<<(n-a+1))-1)*((1ll<<(2*a))-1)/3;
				if (l<=k && k<=r)
				{
					flag=true;
					printf("YES %d\n",a);
					break;
				}
			}
			if (!flag)
			{
				printf("NO\n");
			}
		}
	}
	return 0;
}
