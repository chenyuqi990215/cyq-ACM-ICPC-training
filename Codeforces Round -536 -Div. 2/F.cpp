/*
	原根性质：
		g^(p-1)=1(mod p)
		任意1<=k<p-1，有g^k!=1(mod p)
	推论：
		任意1<=x<y<p-1，有g^x!=g^y(mod p) 
*/
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
struct Matrix
{
	static const int mod=998244352;
	static const int maxm=105;
	ll a[maxm][maxm];
	int row,col;
	Matrix():row(maxm),col(maxm)
	{
		memset(a,0,sizeof(a));
	}
	Matrix(int x,int y):row(x),col(y)
	{
		memset(a,0,sizeof(a));
	}
	ll* operator [](int x)
	{
		return a[x];
	}
	Matrix operator *(Matrix x)
	{
		Matrix tmp(col,x.row);
		for (int i=0;i<row;i++)
			for (int j=0;j<col;j++)
				if (a[i][j])
					for (int k=0;k<x.col;k++)
						if (x[j][k])
							tmp[i][k]=(tmp[i][k]+a[i][j]*x[j][k])%mod;
		return tmp;
	}
	void operator *=(Matrix x)
	{
		*this=*this*x;
	}
	Matrix operator ^(ll x)
	{
		Matrix ret(row,col);
		for (int i=0;i<col;i++)
			ret[i][i]=1;
		Matrix tmp=*this;
		for (;x>0;x>>=1,tmp*=tmp)
			if (x&1) ret*=tmp;
		return ret;
	}
}matrix;
struct BSGS
{
	static const int mod=998244353;
	ll mul_mod(ll a,ll b)
	{
		return a*b%mod;
	}
	ll pow_mod(ll a,ll p)
	{
		if (p==0) return 1;
		ll ans=pow_mod(a,p/2);
		ans=ans*ans%mod;
		if (p&1) ans=ans*a%mod;
		return ans;
	}
	ll inv(ll a)
	{
		return pow_mod(a,mod-2); 
	}
	ll log_mod(ll a,ll b)
	{
		int m=(int)sqrt(mod),v,e=1;
		v=inv(pow_mod(a,m));
		map<int,int> x;
		x[1]=0;
		for (int i=1;i<m;i++)
		{
			e=mul_mod(e,a);
			if (!x.count(e)) x[e]=i;
		}
		for (int i=0;i<m;i++)
		{
			if (x.count(b)) return i*m+x[b];
			b=mul_mod(b,v);
		}
		return -1;
	}
}bsgs;
struct Exgcd
{
	ll gcd(ll a,ll b,ll &x,ll &y)
	{
		if (b==0)
		{
			x=1,y=0;
			return a;
		}
		int q=gcd(b,a%b,y,x);
		y-=a/b*x;
		return q;
	}
	ll gcd(ll a,ll b)
	{
		return b==0?a:gcd(b,a%b);
	}
	//ax+by=c
	bool exgcd(ll a,ll b,ll c,ll &x,ll &y)
	{
		int d=gcd(a,b);
		if (c%d) return false;
		gcd(a,b,x,y);
		while (x<0)
		{
			x+=b/gcd(a,b);
			y+=a/gcd(a,b);
		}
		x*=c/d;
		y*=c/d;
		return true;
	}
}ex;
struct Pow_mod
{
	static const int mod=998244353;
	ll mul_mod(ll a,ll b)
	{
		return a*b%mod;
	}
	ll pow_mod(ll a,ll p)
	{
		if (p==0) return 1;
		ll ans=pow_mod(a,p/2);
		ans=ans*ans%mod;
		if (p&1) ans=ans*a%mod;
		return ans;
	}
}powmod;
int main()
{
	static const int mod=998244353;
	static const int g=3;
	int k,n,m;
	scanf("%d",&k);
	matrix.row=matrix.col=k;
	for (int i=0;i<k;i++)
		scanf("%lld",&matrix.a[0][i]);
	for (int i=1;i<k;i++)
		matrix.a[i][i-1]=1;
	scanf("%d%d",&n,&m);
	matrix=matrix^(n-k);
	ll q=matrix.a[0][0];
	ll p=bsgs.log_mod(g,m);
	ll x,y;
	bool flag=ex.exgcd(q,mod-1,p,x,y);
	if (!flag) printf("-1\n");
	else printf("%lld\n",powmod.pow_mod(3,x));
	return 0;
}
