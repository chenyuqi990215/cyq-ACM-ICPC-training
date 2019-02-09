#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int mod=1e9+7;
struct dp
{
	int f[100000+5];
}g[53],cur[53*6];
int n,m,p;   //m=n/2
int cnt[53],ans[53][53];
void solve(int l,int r,dp &now,int k)
{
	if (l>r) return;
	if (l==r)
	{
		if (cnt[l]<=m) ans[k][l]=ans[l][k]=now.f[m-cnt[l]];
		else ans[k][l]=ans[l][k]=0;
		return;
	}
	int mid=(l+r)/2;
	int nowl=p++;
	int nowr=p++;
	memcpy(cur[nowl].f,now.f,sizeof(now.f));
	for (int i=mid+1;i<=r;i++)
	{
		if (cnt[i]==0) continue;
		for (int j=m;j>=cnt[i];j--)
			cur[nowl].f[j]=(cur[nowl].f[j]+cur[nowl].f[j-cnt[i]])%mod;	
	}	
	solve(l,mid,cur[nowl],k);
	memcpy(cur[nowr].f,now.f,sizeof(now.f));
	for (int i=l;i<=mid;i++)
	{
		if (cnt[i]==0) continue;
		for (int j=m;j>=cnt[i];j--)
			cur[nowr].f[j]=(cur[nowr].f[j]+cur[nowr].f[j-cnt[i]])%mod;	
	}	
	solve(mid+1,r,cur[nowr],k);	
}
int idx(char c)
{
	if (islower(c)) return c-'a'+1;
	else return c-'A'+27;
}
ll pow_mod(ll a,ll p)
{
	if (p==0) return 1;
	ll ret=pow_mod(a,p/2);
	ret=ret*ret%mod;
	if (p%2==1) ret=ret*a%mod;
	return ret;
}
ll inv(ll a)
{
	return pow_mod(a,mod-2);
}
ll fac[100000+5];
void get_fac()
{
	fac[0]=1;
	for (int i=1;i<100000+5;i++)
		fac[i]=(fac[i-1]*i)%mod;
}
char s[100000+5];
int pre[53];
int main()
{
	int q,x,y;
	scanf("%s",s);
	n=strlen(s);
	m=n/2;
	for (int i=0;i<n;i++)
		cnt[idx(s[i])]++;
	g[0].f[0]=1;
	int prei=0;
	for (int i=1;i<=52;i++)
		if (cnt[i])
		{
			pre[i]=prei;
			prei=i;
		}
	int maxi=0;
	for (int i=1;i<=52;i++)
	{
		if (cnt[i]==0) continue;
		maxi=max(maxi,i);
		for (int j=m;j>=0;j--)
		{
			g[i].f[j]=g[pre[i]].f[j];
			if (j>=cnt[i])
				g[i].f[j]=(g[i].f[j]+g[pre[i]].f[j-cnt[i]])%mod;
		}		
		prei=i;
	}
	for (int i=1;i<=52;i++)
	{
		if (cnt[i]==0) continue;
		ans[i][i]=(g[maxi].f[m]*inv(2))%mod;
		p=0;
		int now=p++;
		memset(cur[now].f,0,sizeof(cur[now].f));
		for (int j=cnt[i];j<=m;j++)
			cur[now].f[j]=g[pre[i]].f[j-cnt[i]];
		solve(i+1,52,cur[now],i);
	}
	get_fac();
	ll w=fac[m]*fac[m]%mod;
	for (int i=1;i<=52;i++)
		w=(w*inv(fac[cnt[i]]))%mod;
	for (int i=1;i<=52;i++)
		for (int j=1;j<=52;j++)
		{
			ll tmp=ans[i][j];
			tmp=(tmp+tmp)%mod;
			tmp=(tmp*w)%mod;
			ans[i][j]=(int)tmp;
		}
	scanf("%d",&q);
	for (int i=0;i<q;i++)
	{
		scanf("%d%d",&x,&y);
		printf("%d\n",ans[idx(s[x-1])][idx(s[y-1])]);
	}
	return 0;
}
