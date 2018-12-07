#include <bits/stdc++.h>
using namespace std;
const int maxn=1e5+5;
const int maxm=1e6+5;
const int mod=1e9+7;
vector <int> fac[maxm];
int a[maxn];
long long dp[maxm];
void get_fac()
{
	for (int i=1;i<maxm;i++)
		for (int j=i;j<maxm;j+=i)
			fac[j].push_back(i);
}
int main()
{
	int n; 
	get_fac();
	scanf("%d",&n);
	for (int i=0;i<n;i++)
		scanf("%d",&a[i]);
	memset(dp,0,sizeof(dp));
	dp[0]=1;
	for (int i=0;i<n;i++)
	{
		int k=a[i];
		for (int j=fac[k].size()-1;j>=0;j--)
			dp[fac[k][j]]=(dp[fac[k][j]]+dp[fac[k][j]-1])%mod;
	}
	long long ans=0;
	for (int i=1;i<maxn;i++)
		ans=(ans+dp[i])%mod;
	printf("%d\n",(int)ans);
	return 0;
}
