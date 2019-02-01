#include <bits/stdc++.h>
int a[105];
bool mark[105];
using namespace std;
int main()
{
	int n,k;
	scanf("%d%d",&n,&k);
	for (int i=1;i<=n;i++)
		scanf("%d",&a[i]);
	int ans=0;
	for (int b=1;b<=n;b++)
	{
		int t=b,u=b;
		memset(mark,0,sizeof(mark));
		while (t>=1)
		{
			mark[t]=1;
			t-=k;
		}
		while (u<=n)
		{
			mark[u]=1;
			u+=k;
		}
		int e=0,s=0;
		for (int i=1;i<=n;i++)
		{
			if (a[i]==1 && mark[i]==0) e++;
			if (a[i]==-1 && mark[i]==0) s++;
		}
		ans=max(ans,abs(e-s));
	}
	printf("%d\n",ans);
	return 0;
}
