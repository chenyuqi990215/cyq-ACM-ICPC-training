#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int maxn=4*1000000+5;
ll c[2][maxn];
vector <int> g[100000+5];
int t[100000+5],p[100000+5],n;
ll x[100000+5],candy[100000+5],dp[100000+5],l[100000+5],T;
int lowbit(int x)
{
	return x&(-x);
}
ll sum(int x,int p)
{
	ll ret=0;
	while (x>0)
	{
		ret+=c[p][x];
		x-=lowbit(x);
	}
	return ret;
}
void update(int x,int p,ll d)
{
	while (x<maxn)
	{
		c[p][x]+=d;
		x+=lowbit(x);
	}
}
void dfs(int u,ll k)
{
	k+=2ll*l[u];
	ll left=T-k;
	update(t[u],0,1ll*t[u]*x[u]);
	update(t[u],1,x[u]);
	ll l=1,r=1000000;
	while (l<=r)
	{
		ll mid=(l+r)/2;
		if (sum(mid,0)>left) r=mid-1;
		else l=mid+1;
	}
	candy[u]=sum(r,1);
	left-=sum(r,0);
	candy[u]+=max(0ll,min(left/l,sum(l,1)-sum(r,1)));
	dp[u]=candy[u];
	ll maxc=-1,semimaxc=-1;
	for (int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		dfs(v,k);
		if (dp[v]>maxc)
		{
			semimaxc=maxc;
			maxc=dp[v];
		}
		else if (dp[v]>semimaxc)
			semimaxc=dp[v];
	}
	dp[u]=max(dp[u],semimaxc);
	update(t[u],0,-1ll*t[u]*x[u]);
	update(t[u],1,-x[u]);
}
int main()
{
	scanf("%d%lld",&n,&T);
	for (int i=1;i<=n;i++)
		scanf("%lld",&x[i]);
	for (int i=1;i<=n;i++)
		scanf("%d",&t[i]);
	for (int i=2;i<=n;i++)
	{
		scanf("%d%lld",&p[i],&l[i]);
		g[p[i]].push_back(i);
	}
	dfs(1,0);
	ll ans=dp[1];
	for (int i=0;i<g[1].size();i++)
	{
		int v=g[1][i];
		ans=max(ans,dp[v]);
	}
	printf("%lld\n",ans);
	return 0;
}
