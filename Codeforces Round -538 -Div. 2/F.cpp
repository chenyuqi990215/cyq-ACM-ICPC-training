#include <bits/stdc++.h>
#define lc (2*o)
#define rc (2*o+1)
#define mp make_pair
using namespace std;
const int mod=1e9+7;
const int maxn=2e6+5;
typedef pair<int,int> pii;
typedef long long ll;
struct node
{
	int op,l,r,x;
	ll ans;
}qnode[200000+5],anode[400000+5];
ll pow_mod(ll a,ll p)
{
	if (p==0) return 1;
	ll ret=pow_mod(a,p/2);
	ret=ret*ret%mod;
	if (p%2==1) ret=ret*a%mod;
	return ret;
} 
ll _div[305];
vector <int> is_prime;
bool prime[305];
ll inv[65];
void get_div()
{
	is_prime.push_back(2);
	memset(prime,true,sizeof(prime));
	for (int i=3;i<=300;i+=2)
		if (prime[i])
		{
			is_prime.push_back(i);
			for (int j=2*i;j<=300;j+=i)
				prime[j]=false;
		}
	memset(_div,0,sizeof(_div));
	for (int i=1;i<=300;i++)
	{
		for (int j=0;j<is_prime.size();j++)
		{
			if (i%is_prime[j]==0)
				_div[i]|=1ll<<j;
		}
	}
	for (int i=0;i<is_prime.size();i++)
		inv[i]=pow_mod(is_prime[i],mod-2);
}
struct segmul
{
	ll sumv[maxn],mulv[maxn],v1;
	ll bitv[maxn],sbitv[maxn],v2;
	int n,Y1,Y2;
	void build(int o,int l,int r)
	{
		if (l==r)
		{
			mulv[o]=anode[l].x;
			sumv[o]=anode[l].x;
			bitv[o]=_div[anode[l].x];
			sbitv[o]=_div[anode[l].x];
		}
		else
		{
			int mc=(l+r)>>1;
			build(lc,l,mc);
			build(rc,mc+1,r);
			sumv[o]=sumv[lc]*sumv[rc]%mod;
			bitv[o]=bitv[lc]|bitv[rc];
			mulv[o]=1,sbitv[o]=0;
		}
	}
	void build()
	{
		build(1,1,n);
	}
	void pushup(int o)
	{
		sumv[o]=sumv[lc]*sumv[rc]%mod;
		bitv[o]=bitv[lc]|bitv[rc];
	}
	void pushdown(int o,int m)
	{
		mulv[lc]=(mulv[lc]*mulv[o])%mod;
		mulv[rc]=(mulv[rc]*mulv[o])%mod;
		sumv[lc]=(sumv[lc]*pow_mod(mulv[o],(m-(m>>1))))%mod;
		sumv[rc]=(sumv[rc]*pow_mod(mulv[o],(m>>1)))%mod;
		mulv[o]=1;
		bitv[lc]|=sbitv[o];
		bitv[rc]|=sbitv[o];
		sbitv[lc]|=sbitv[o];
		sbitv[rc]|=sbitv[o];
		sbitv[o]=0;
	}
	void _update(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
		{
			sumv[o]=(sumv[o]*pow_mod(v1,(r-l+1)))%mod;
			mulv[o]=(mulv[o]*v1)%mod;
			bitv[o]|=v2;
			sbitv[o]|=v2;
		}
		else
		{
			pushdown(o,r-l+1);
			int mc=(l+r)>>1;
			if (Y1<=mc) _update(lc,l,mc);
			if (Y2>mc) _update(rc,mc+1,r);
			pushup(o);
		}	
	}
	pair<ll,ll> _query(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
			return mp(sumv[o],bitv[o]);
		else
		{
			pushdown(o,r-l+1);
			int mc=(l+r)>>1;
			ll ansmul=1,ansbit=0;
			if (Y1<=mc) 
			{
				pair<ll,ll> tmp=_query(lc,l,mc);
				ansmul=(ansmul*tmp.first)%mod;
				ansbit|=tmp.second;
			}
			if (Y2>mc)
			{
				pair<ll,ll> tmp=_query(rc,mc+1,r);
				ansmul=(ansmul*tmp.first)%mod;
				ansbit|=tmp.second;
			}
			return mp(ansmul,ansbit);
		}
	}
	void update(int Y1,int Y2,ll v1,ll v2)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		this->v1=v1;
		this->v2=v2;
		_update(1,1,n); 
	}
	pair<ll,ll> query(int Y1,int Y2)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		return _query(1,1,n);
	}
}tr;
int main()
{
	int n,q;
	scanf("%d%d",&n,&q);
	get_div();
	tr.n=n;
	for (int i=1;i<=n;i++)
	{
		scanf("%d",&anode[i].x);
		anode[i].op=1;
		anode[i].l=anode[i].r=i;
	}
	char s[20];
	for (int i=0;i<q;i++)
	{
		scanf("%s",s);
		qnode[i].ans=1;
		if (strcmp(s,"MULTIPLY")==0)
		{
			qnode[i].op=1;
			scanf("%d%d%d",&qnode[i].l,&qnode[i].r,&qnode[i].x);
		}
		else
		{
			qnode[i].op=2;
			scanf("%d%d",&qnode[i].l,&qnode[i].r);
		}
	}
	tr.build();
	for (int i=0;i<q;i++)
	{
		if (qnode[i].op==1)
			tr.update(qnode[i].l,qnode[i].r,qnode[i].x,_div[qnode[i].x]);
		else
		{
			pair<ll,ll> cur=tr.query(qnode[i].l,qnode[i].r);
			qnode[i].ans=cur.first;
			ll tmp=cur.second;
			for (int j=0;j<is_prime.size();j++)
			{
				if (tmp&(1ll<<j))
				{
					qnode[i].ans=(qnode[i].ans*inv[j])%mod;
					qnode[i].ans=(qnode[i].ans*(is_prime[j]-1))%mod;
				}	
			}
		}
	}
	for (int i=0;i<q;i++)
	{
		if (qnode[i].op==2)
			printf("%lld\n",qnode[i].ans);
	}
	return 0;
}
