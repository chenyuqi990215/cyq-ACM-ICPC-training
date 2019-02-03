#include <bits/stdc++.h>
using namespace std;
vector <int> g[100000+5];
bool done[100000+5];
int main()
{
	int n,m,u,v;
	scanf("%d%d",&n,&m);
	for (int i=0;i<m;i++)
	{
		scanf("%d%d",&u,&v);
		g[u].push_back(v);
		g[v].push_back(u);
	}
	priority_queue <int,vector <int>, greater<int> >pq;
	vector <int> ans;
	pq.push(1);
	while (!pq.empty())
	{
		int x=pq.top();
		pq.pop();
		if (done[x]) continue;
		ans.push_back(x);
		done[x]=true;
		for (int i=0;i<g[x].size();i++)
			pq.push(g[x][i]);
	}
	printf("%d",ans[0]);
	for (int i=1;i<ans.size();i++)
		printf(" %d",ans[i]);
	printf("\n");
	return 0;
}
