#include <bits/stdc++.h>
#define mp make_pair
using namespace std;
const int maxn=1e6+5;
typedef long long ll;
typedef pair<ll,int> pli; 
vector <int> is_prime;
bool prime[maxn]; 
vector <pli> divi;
vector <ll> cnt;
int main()
{
	ll n,b;
	scanf("%lld%lld",&n,&b);
	memset(prime,true,sizeof(prime));
	is_prime.push_back(2);
	for (int i=4;i<maxn;i+=2)
		prime[i]=false;
	for (int i=3;i<maxn;i+=2)
	{
		if (prime[i])
		{
			is_prime.push_back(i);
			for (int j=2*i;j<maxn;j+=i)
				prime[j]=false;
		}
	}
	for (int i=0;i<is_prime.size();i++)
	{
		if (b%is_prime[i]==0)
		{
			int cnt=0;
			while (b%is_prime[i]==0)
				b/=is_prime[i],cnt++;
			divi.push_back(mp(is_prime[i],cnt));
		}
	}
	if (b!=1)
		divi.push_back(mp(b,1));
	ll ans=-1;
	for (int i=0;i<divi.size();i++)
	{
		cnt.clear();
		ll tmp=n;
		ll cur=divi[i].first;
		while (tmp)
			cnt.push_back(tmp/cur),tmp/=cur;
		for (int j=cnt.size()-1;j>=0;j--)
			for (int k=j+1;k<cnt.size();k++)
				cnt[j]-=cnt[k];
		ll sum=0;
		for (int j=0;j<cnt.size();j++)
			sum+=1ll*cnt[j]*(j+1);
		sum/=divi[i].second;
		if (ans==-1) ans=sum;
		else ans=min(ans,sum);
	}
	printf("%lld\n",ans);
	return 0;
}
