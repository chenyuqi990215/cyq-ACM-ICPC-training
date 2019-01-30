#include <bits/stdc++.h>
using namespace std;
const int maxn=300000+5;
const int INF=2e9;
struct node
{
	int l,r;
	int minv;
}T[maxn*30];
int rt[maxn];
int tot=0;
void build(int &o,int l,int r)
{
	o=++tot;
	T[o].minv=INF;
	if (l==r) return;
	int mid=(l+r)/2;
	build(T[o].l,l,mid);
	build(T[o].r,mid+1,r);
}
void update(int l,int r,int &now,int last,int p,int v)
{
	T[++tot]=T[last];
	now=tot;
	if (l==r)
	{
		T[now].minv=min(T[now].minv,v);
		return;
	}
	int mid=(l+r)/2;
	if (p<=mid)
		update(l,mid,T[now].l,T[last].l,p,v);
	else update(mid+1,r,T[now].r,T[last].r,p,v);
	T[now].minv=max(T[T[now].l].minv,T[T[now].r].minv);
}
int query(int o,int l,int r,int y1,int y2)
{
	if (y1<=l && r<=y2) return T[o].minv;
	int ans=0;
	int mid=(l+r)/2;
	if (y1<=mid) ans=max(ans,query(T[o].l,l,mid,y1,y2));
	if (y2>mid) ans=max(ans,query(T[o].r,mid+1,r,y1,y2));
	return ans;
}
struct data
{
	int l,r,p;
	bool operator < (const data &rhs) const
	{
		return l>rhs.l;
	}
}seg[maxn];
int main()
{
	int n,m,k,a,b,x,y;
	scanf("%d%d%d",&n,&m,&k);
	for (int i=1;i<=k;i++)
		scanf("%d%d%d",&seg[i].l,&seg[i].r,&seg[i].p);
	sort(seg+1,seg+k+1);
	build(rt[0],1,n);
	for (int i=1;i<=k;i++)
		update(1,n,rt[i],rt[i-1],seg[i].p,seg[i].r);
	for (int i=0;i<m;i++)
	{
		scanf("%d%d%d%d",&a,&b,&x,&y);
		int t=(upper_bound(seg+1,seg+k+1,data{x,-1,-1})-seg)-1;
		int maxv=query(rt[t],1,n,a,b);
		if (maxv>y) printf("no\n");
		else printf("yes\n");
		fflush(stdout);
	}
	return 0;
}
