#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
ll cnt[3],dp[200005][3];
const ll mod=1e9+7;
void cal(int x,int sign)
{
	if (x==0) return;
	cnt[0]+=sign*(x/3);
	cnt[1]+=sign*((x+2)/3);
	cnt[2]+=sign*((x+1)/3);
}
int main()
{
	int n,l,r;
	scanf("%d%d%d",&n,&l,&r);
	memset(dp,0,sizeof(dp));
	dp[0][0]=1;
	cal(l-1,-1);
	cal(r,1);	
	for (int i=1;i<=n;i++)
		for (int j=0;j<3;j++)
			for (int k=0;k<3;k++)
				dp[i][(j+k)%3]=(dp[i][(j+k)%3]+dp[i-1][j]*cnt[k])%mod;
	printf("%lld\n",dp[n][0]);
	return 0;
}
