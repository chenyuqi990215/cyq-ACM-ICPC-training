#include <bits/stdc++.h>
#define mp make_pair
#define lc 2*o
#define rc 2*o+1
#define mc l+(r-l)/2
using namespace std;
const int maxn=100000+5;
typedef long long ll;
typedef pair<int,int> pii;
struct Segtree
{
	pii maxv[4*maxn],v;
	int p,ql,qr,n;
	void init(int n)
	{
		this->n=n;
		for (int i=0;i<=4*n;i++)
			maxv[i]=mp(-1,-1);
	}
	void update(int o,int l,int r)
	{
		if (l==r)
			maxv[o]=max(maxv[o],v);
		else
		{
			if (p<=mc) update(lc,l,mc);
			else update(rc,mc+1,r);
			maxv[o]=max(maxv[lc],maxv[rc]);
		}
	}
	pii quiry(int o,int l,int r)
	{
		pii ans=mp(-1,-1);
		if (ql<=l && r<=qr) return maxv[o];
		if (ql<=mc) ans=max(ans,quiry(lc,l,mc));
		if (qr>mc) ans=max(ans,quiry(rc,mc+1,r));
		return ans;
	}
	void update(int p,pii v)
	{
		this->p=p;
		this->v=v;
		update(1,1,n);
	}
	pii quiry(int ql,int qr)
	{
		this->ql=ql;
		this->qr=qr;
		return quiry(1,1,n);
	}
}tr;
struct node
{
	int s,t,d,w;
	bool operator < (const node &rhs) const
	{
		return s<rhs.s;
	} 
}data[maxn];
ll dp[maxn][205];
pii ch[maxn];
int main()
{
	int n,m,k;
	scanf("%d%d%d",&n,&m,&k);
	tr.init(n);
	for (int i=0;i<k;i++)
		scanf("%d%d%d%d",&data[i].s,&data[i].t,&data[i].d,&data[i].w);
	sort(data,data+k);
	int p=0;
	for (int i=1;i<=n;i++)
	{
		while (p<k && data[p].s<=i) 
		{
			tr.update(data[p].t,mp(data[p].w,data[p].d));
			p++;
		}
		ch[i]=tr.quiry(i,n);
	}
	memset(dp,0,sizeof(dp));
	if (ch[n]!=mp(-1,-1))
		dp[n][0]=ch[n].first;
	for (int i=n-1;i>=1;i--)
		for (int j=0;j<=m;j++)
		{
			dp[i][j]=1e15;
			if (ch[i].first!=-1)
			{
				dp[i][j]=min(dp[i][j],dp[ch[i].second+1][j]+ch[i].first);  //不打扰 
				if (j>0) dp[i][j]=min(dp[i][j],dp[i+1][j-1]);   //打扰 
			}
			else
				dp[i][j]=dp[i+1][j];   //不打扰 
		}
	ll ans=1e15;
	for (int i=0;i<=m;i++)
		ans=min(ans,dp[1][i]);
	printf("%lld\n",ans);
	return 0;
}
