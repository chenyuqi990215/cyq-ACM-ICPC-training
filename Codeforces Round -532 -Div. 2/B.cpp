#include <bits/stdc++.h>
using namespace std;
int minv[4*100000+5];
void update(int o,int l,int r,int p,int v)
{
	int lc=2*o,rc=2*o+1;
	if (l==r)
	{
		minv[o]+=v;
		return;
	}
	int m=l+(r-l)/2;
	if (p<=m) update(lc,l,m,p,v);
	else update(rc,m+1,r,p,v);
	minv[o]=min(minv[lc],minv[rc]);
}
int main()
{
	int n,m,x;
	scanf("%d%d",&n,&m);
	int k=1;
	for (int i=0;i<m;i++)
	{
		scanf("%d",&x);
		update(1,1,n,x,1);
		if (minv[1]==k)
		{
			printf("1");
			k++;
		}
		else
			printf("0");
	}
	return 0;
}
