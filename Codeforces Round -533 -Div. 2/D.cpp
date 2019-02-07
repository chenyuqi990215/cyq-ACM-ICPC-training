#include <bits/stdc++.h>
#define mp make_pair
using namespace std;
struct node
{
	int x,y,d;
};
int dirx[]={1,-1,0,0};
int diry[]={0,0,1,-1};
int n,m,p,s[10],a[1005][1005],cnt[10];
char t[1005][1005];
void bfs()
{
	queue <node> q[10][2];
	int r=0;
	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			if (a[i][j]>0)
				q[a[i][j]][r].push(node{i,j,0});
	while (true)
	{
		bool finish=true;
		for (int i=1;i<=9;i++)
			if (!q[i][r].empty())
				finish=false;
		if (finish) break;
		for (int i=1;i<=p;i++)
		{
			while (!q[i][r^1].empty())
				q[i][r^1].pop();
			while (!q[i][r].empty())
			{
				node p=q[i][r].front();
				q[i][r].pop();
				int x=p.x,y=p.y,d=p.d;
				if (d==s[i])
				{
					q[i][r^1].push(node{x,y,0});
					continue;
				} 
				for (int k=0;k<4;k++)
				{
					int nx=x+dirx[k],ny=y+diry[k];
					if (nx<0 || nx>=n || ny<0 || ny>=m) continue;
					if (a[nx][ny]!=0) continue;
					a[nx][ny]=i;
					q[i][r].push(node{nx,ny,d+1});
				}
			}
		}
		r^=1;
	}
}
int main()
{
	scanf("%d%d%d",&n,&m,&p);
	for (int i=1;i<=p;i++)
		scanf("%d",&s[i]);
	for (int i=0;i<n;i++)
		scanf("%s",t[i]);
	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
		{
			if (t[i][j]=='#') a[i][j]=-1;
			else if (t[i][j]=='.') a[i][j]=0;
			else a[i][j]=t[i][j]-'0';
		}
	bfs();
	for (int i=0;i<n;i++)
		for (int j=0;j<m;j++)
			if (a[i][j]>0)
				cnt[a[i][j]]++;
	for (int i=1;i<=p;i++)
		printf("%d%c",cnt[i],i==p?'\n':' ');
	return 0;
}
