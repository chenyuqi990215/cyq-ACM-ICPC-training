#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int a[100000+5],n,k,A,B;
ll solve(int l,int r)
{
	int na=upper_bound(a,a+k,r)-lower_bound(a,a+k,l);
	if (na==0) return A;
	ll ret=1ll*B*na*(r-l+1);
	if (l==r) return ret;
	int mid=(l+r)/2;
	ret=min(ret,solve(l,mid)+solve(mid+1,r));
	return ret;
}
int main()
{
	scanf("%d%d%d%d",&n,&k,&A,&B);
	for (int i=0;i<k;i++)
		scanf("%d",&a[i]);
	sort(a,a+k);
	printf("%lld\n",solve(1,1<<n));
	return 0;
}
