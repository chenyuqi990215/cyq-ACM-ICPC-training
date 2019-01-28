#include <bits/stdc++.h>
using namespace std;
int main()
{
	int w,h,u1,d1,u2,d2;
	scanf("%d%d%d%d%d%d",&w,&h,&u1,&d1,&u2,&d2);
	for (int i=h;i>=0;i--)
	{
		w+=i;
		if (i==d1)
			w=max(0,w-u1);
		else if (i==d2)
			w=max(0,w-u2);
	}
	printf("%d\n",w);
	return 0;
} 
