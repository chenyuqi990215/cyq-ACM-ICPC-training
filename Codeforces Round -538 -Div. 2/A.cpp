#include <bits/stdc++.h>
using namespace std;
int main()
{
	int a,b,c,x,y,z;
	scanf("%d%d%d%d%d%d",&x,&y,&z,&a,&b,&c);
	bool flag=(a>=x && a+b>=x+y && a+b+c>=x+y+z);
	printf("%s\n",flag?"YES":"NO");
	return 0; 
}
