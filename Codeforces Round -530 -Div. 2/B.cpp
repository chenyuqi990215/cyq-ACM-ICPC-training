#include <bits/stdc++.h>
using namespace std;
int main()
{
	int n;
	scanf("%d",&n);
	for (int i=1,k;;i++)
	{
		if (i%2==1)
			k=((i-1)/2)*((i+1)/2);
		else k=(i/2)*(i/2);
		if (k>=n)
		{
			printf("%d\n",i);
			break;
		}
	}
	return 0;
}
