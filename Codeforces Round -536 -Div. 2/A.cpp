#include <bits/stdc++.h>
using namespace std;
int n;
char s[505][505];
int main()
{
	scanf("%d",&n);
	for (int i=0;i<n;i++)
		scanf("%s",s[i]);
	int cnt=0;
	for (int i=1;i<n-1;i++)
		for (int j=1;j<n-1;j++)
		{
			int temp=0;
			temp+=(s[i][j]=='X');
			temp+=(s[i-1][j-1]=='X');
			temp+=(s[i-1][j+1]=='X');
			temp+=(s[i+1][j-1]=='X');
			temp+=(s[i+1][j+1]=='X');
			if (temp==5) cnt++;
		}
	printf("%d\n",cnt);
	return 0;
} 
