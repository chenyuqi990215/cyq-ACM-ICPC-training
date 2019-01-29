#include <bits/stdc++.h>
using namespace std;
int cal(int x)
{
	if (x%2==0) return x/2;
	return (x-1)/2-x;
}
int main()
{
	int t,l,r;
	scanf("%d",&t);
	while (t--)
	{
		scanf("%d%d",&l,&r);
		printf("%d\n",cal(r)-cal(l-1));
	}
	return 0;
}
