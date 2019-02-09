#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int mod=1e9+7;
vector <int> g[100000+5];
const int LN=20;
int par[100000+5][20],lev[100000+5];
int st[100000+5],en[100000+5],dfs_clock=0;
void dfs(int u,int p)
{
	st[u]=++dfs_clock;
    if (p!=-1)
    {
        lev[u]=lev[p]+1;
        par[u][0]=p;
        for (int j=1;j<LN;j++)
            if (par[u][j-1]!=-1)
                par[u][j]=par[par[u][j-1]][j-1];
    }
    for (int i=0;i<g[u].size();i++)
    {
    	int v=g[u][i];
    	if (v==p) continue;
    	dfs(v,u);
	}
	en[u]=++dfs_clock;
}
int lca(int u,int v)
{
    if (lev[u]<lev[v])
        swap(u,v);
    for (int j=LN-1;j>=0;j--)
        if (par[u][j]+1 && lev[par[u][j]]>=lev[v])
            u=par[u][j];
    if (u==v)
        return u;
    for (int j=LN-1;j>=0;j--)
    {
        if (par[u][j]+1 && par[u][j]!=par[v][j])
        {
            u=par[u][j];
            v=par[v][j];
        }
    }
    return par[u][0];
}
int h[100000+5];
int lowbit(int x)
{
	return x&(-x);
}
int c[800000+5];
int sum(int x)
{
	int ret=0;
	while (x>0)
	{
		ret+=c[x];
		x-=lowbit(x);
	}
	return ret;
}
void update(int x,int v)
{
	while (x<800000+5)
	{
		c[x]+=v;
		x+=lowbit(x);
	}
}
int a[100000+5];
bool mark[100000+5];
ll dp[100000+5];
int main()
{
	int n,q,k,m,r,u,v;
	scanf("%d%d",&n,&q);
	for (int i=1;i<n;i++)
	{
		scanf("%d%d",&u,&v);
		g[u].push_back(v);
		g[v].push_back(u);
	}
	memset(par,-1,sizeof(par));
	dfs(1,-1);
	while (q--)
	{
		scanf("%d%d%d",&k,&m,&r);
		for (int i=0;i<k;i++)
		{
			scanf("%d",&a[i]);
			mark[a[i]]=true;
		}
		for (int i=0;i<k;i++)
		{
			update(st[a[i]],1);
			update(en[a[i]]+1,-1);
		}
		for (int i=0;i<k;i++)
		{
			int _lca=lca(r,a[i]);
			h[i]=sum(st[a[i]])+sum(st[r])-2*sum(st[_lca])+(int)mark[_lca]-1;
		}
		for (int i=0;i<k;i++)
		{
			update(st[a[i]],-1);
			update(en[a[i]]+1,1);
			mark[a[i]]=false;
		}
		sort(h,h+k);
		for (int i=0;i<=m;i++)
			dp[i]=0;
		dp[0]=1;
		for (int i=0;i<k;i++)
			for (int j=m;j>=0;j--)
			{
				if (j<=h[i]) dp[j]=0;
				else dp[j]=(dp[j]*(j-h[i])+dp[j-1])%mod;
			}
		ll ans=0;
		for (int i=1;i<=m;i++)
			ans=(ans+dp[i])%mod;
		printf("%lld\n",ans);
	}
	return 0;
}
