#include <bits/stdc++.h>
using namespace std;
int n,k;
char s[200000+5];
int cal(char c)
{
	int ret=0,p=0;
	while (p<n)
	{
		if (s[p]==c)
		{
			int t=0;
			for (;t<k && p<n;p++,t++)
				if (s[p]!=c) break;
			if (t==k) ret++;
			p--;
		}
		p++;
	}
	return ret;
}
int main()
{
	scanf("%d%d",&n,&k);
	scanf("%s",s);
	int ans=0;
	for (char c='a';c<='z';c++)
		ans=max(ans,cal(c));
	printf("%d\n",ans);
	return 0;
}
