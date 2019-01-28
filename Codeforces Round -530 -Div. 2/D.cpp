#include <bits/stdc++.h>
using namespace std;
int s[100000+5],sum[100000+5],a[100000+5];
const int INF=2e9;
vector <int> g[100000+5];
bool flag=true;
void dfs(int u,int fa,int dep)
{
	if (dep==1)
	{
		if (s[u]==-1)
			flag=false;
		a[u]=s[u];
		sum[u]=a[u];
	}
	else if (dep%2==1)
	{
		if (s[u]==-1)
			flag=false;
		if (s[u]<sum[fa])
			flag=false;
		else 
		{
			a[u]=s[u]-sum[fa];
			sum[u]=sum[fa]+a[u];
		}
	}
	else
	{
		if (s[u]!=-1)
			flag=false;
		int minx=INF;
		for (int i=0;i<g[u].size();i++)
		{
			int v=g[u][i];
			minx=min(minx,s[v]);
			if (minx<s[fa])
				flag=false;
			a[u]=minx-s[fa];
			sum[u]=sum[fa]+a[u];
		}
	}		
	for (int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		dfs(v,u,dep+1);
	}
}
int main()
{
	int n,x;
	scanf("%d",&n);
	for (int i=2;i<=n;i++)
	{
		scanf("%d",&x);
		g[x].push_back(i);
	} 
	for (int i=1;i<=n;i++)
		scanf("%d",&s[i]);
	long long ans=0;
	dfs(1,-1,1);
	if (!flag) printf("-1\n");
	else
	{
		for (int i=1;i<=n;i++)
			ans+=a[i];
		printf("%lld\n",ans);
	}
	return 0;
}
