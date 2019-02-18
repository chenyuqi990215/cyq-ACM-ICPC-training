#include <bits/stdc++.h>
#define mp make_pair
using namespace std;
typedef pair<int,int> pii;
typedef long long ll;
pii a[200000+5];
int b[200000+5];
int main()
{
	int n,m,k,x;
	scanf("%d%d%d",&n,&m,&k);
	for (int i=1;i<=n;i++)
	{
		scanf("%d",&x);
		a[i]=mp(x,i);
	}
	sort(a+1,a+n+1);
	ll sum=0;
	for (int i=n,j=1;j<=m*k;i--,j++)
	{
		sum+=a[i].first;
		b[j]=a[i].second;
	}
	printf("%lld\n",sum);
	sort(b+1,b+m*k+1);
	printf("%d",b[m]);
	for (int i=2*m;i<m*k;i+=m)
		printf(" %d",b[i]);
	printf("\n");
	return 0;
}
