#include <bits/stdc++.h>
using namespace std;
struct Edge
{
	int u,v,c;
}edge[100000+5];
int n,m;
vector <int> g[100000+5];
int ord[100000+5];
int in[100000+5];
bool check(int x)
{
	memset(in,0,sizeof(in));
	for (int i=1;i<=m;i++)
	{
		if (edge[i].c>x)
			in[edge[i].v]++;
	}
	queue <int> q;
	int cnt=0;
	for (int i=1;i<=n;i++)
		if (in[i]==0)
		{
			q.push(i);
			ord[i]=cnt;
			cnt++;
		}
	while (!q.empty())
	{
		int u=q.front();
		q.pop();
		for (int i=0;i<g[u].size();i++)
		{
			if (edge[g[u][i]].c>x)
			{
				in[edge[g[u][i]].v]--;
				if (in[edge[g[u][i]].v]==0)
				{
					q.push(edge[g[u][i]].v);
					ord[edge[g[u][i]].v]=cnt;
					cnt++;
				} 
			}
		}
	}
	return cnt==n;
} 
int main()
{
	scanf("%d%d",&n,&m);
	for (int i=1;i<=m;i++)
	{
		scanf("%d%d%d",&edge[i].u,&edge[i].v,&edge[i].c);
		g[edge[i].u].push_back(i);
	}
	int l=0,r=1e9;
	while (l<=r)
	{
		int mid=(l+r)/2;
		if (check(mid)) r=mid-1;
		else l=mid+1;
	}
	check(l);
	vector <int> s;
	for (int i=1;i<=m;i++)
	{
		if (ord[edge[i].u]>ord[edge[i].v])
			s.push_back(i);
	}
	printf("%d %d\n",l,s.size());
	for (int i=0;i<s.size();i++)
	{
		printf("%d%c",s[i],(i==s.size()-1)?'\n':' ');
	}
	return 0;
}
