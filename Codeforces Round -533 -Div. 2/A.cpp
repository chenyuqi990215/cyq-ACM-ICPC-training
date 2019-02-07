#include <bits/stdc++.h>
using namespace std;
int a[1005];
int main()
{
	int n;
	int best=-1;
	scanf("%d",&n);
	for (int i=0;i<n;i++)
		scanf("%d",&a[i]);
	int ans=100*n;
	for (int t=1;t<=101;t++)
	{
		int tmp=0;
		for (int i=0;i<n;i++)
		{
			int cur=0;
			if (a[i]>t+1) cur=a[i]-t-1;
			if (a[i]<t-1) cur=t-1-a[i];
			tmp+=cur;
		}
		if (tmp<ans)
		{
			ans=tmp;
			best=t;
		}
	}
	printf("%d %d\n",best,ans);
	return 0;
}
