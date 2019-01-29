#include <bits/stdc++.h>
using namespace std;
int main()
{
	int n,k;
	scanf("%d%d",&n,&k);
	int ans=0;
	ans+=(2*n-1)/k+1;
	ans+=(5*n-1)/k+1;
	ans+=(8*n-1)/k+1;
	printf("%d\n",ans);
	return 0;
} 
