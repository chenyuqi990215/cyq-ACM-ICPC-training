#include <bits/stdc++.h>
using namespace std;
int n,dp[5005][5005][2],a[5005];
int solve(int l,int r,int s)
{
	if (l==r) return dp[l][r][s]=0;
	if (dp[l][r][s]>=0) return dp[l][r][s];
	int &ret=dp[l][r][s];
	ret=6e3;
	if (s==0)
	{
		if (l<n)
			ret=min(ret,solve(l+1,r,0)+(a[l+1]!=a[l]));
		if (l<n)
			ret=min(ret,solve(l+1,r,1)+(a[l]!=a[r]));
		if (r>1 && a[l]==a[r])
			ret=min(ret,solve(l,r-1,0)); 
		if (r>1 && a[r-1]==a[r])
			ret=min(ret,solve(l,r-1,1)+(a[l]!=a[r-1]));
	}
	else
	{
		if (r>1)
			ret=min(ret,solve(l,r-1,1)+(a[r-1]!=a[r]));
		if (r>1)
			ret=min(ret,solve(l,r-1,0)+(a[l]!=a[r]));
		if (l<n && a[r]==a[l])
			ret=min(ret,solve(l+1,r,1));
		if (l<n && a[l+1]==a[l])
			ret=min(ret,solve(l+1,r,0)+(a[l+1]!=a[r]));
	}
	return ret;
}
int main()
{
	memset(dp,-1,sizeof(dp));
	scanf("%d",&n);
	for (int i=1;i<=n;i++)
		scanf("%d",&a[i]);
	printf("%d\n",min(solve(1,n,0),solve(1,n,1)));
	return 0;
}
