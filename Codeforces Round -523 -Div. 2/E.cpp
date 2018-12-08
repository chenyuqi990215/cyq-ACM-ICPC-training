#include <bits/stdc++.h>
using namespace std;
const int maxn=1e3+50;
const int INF=1e9+7;
struct Edge
{
	int from,to,cap,flow,cost;
	Edge(int from,int to,int cap,int flow,int cost):from(from),to(to),cap(cap),flow(flow),cost(cost){}
}; 
struct MCMF
{
	vector <Edge> edge;
	vector <int> g[maxn];
	int a[maxn],p[maxn],c[maxn],inq[maxn],n;
	void init(int n)
	{
		this->n=n;
		edge.clear();
		for (int i=0;i<=n;i++)
			g[i].clear();
	}
	void addedge(int from,int to,int cap,int cost)
	{
		edge.push_back(Edge(from,to,cap,0,cost));
		edge.push_back(Edge(to,from,0,0,-cost));
		int m=edge.size();
		g[from].push_back(m-2);
		g[to].push_back(m-1);
	}
	bool bellmanford(int s,int t,int limit_flow,int &flow,long long &cost)
	{
		memset(a,0,sizeof(a));
		memset(c,126,sizeof(c));
		memset(inq,0,sizeof(inq));
		queue <int> q;
		a[s]=INF;
		c[s]=0;
		inq[s]=1;
		q.push(s);
		while (!q.empty())
		{
			int x=q.front();
			inq[x]=0;
			q.pop();
			for (int i=0;i<g[x].size();i++)
			{
				Edge &e=edge[g[x][i]];
				if (e.cap>e.flow && c[e.to]>c[x]+e.cost)
				{
					c[e.to]=c[x]+e.cost;
					p[e.to]=g[x][i];
					a[e.to]=min(a[x],e.cap-e.flow);
					if (!inq[e.to])
					{
						q.push(e.to);
						inq[e.to]=1;
					}
				}
			}
		} 
		if (c[t]>=INF) return false;
		if (a[t]+flow>limit_flow) a[t]=limit_flow-flow;
		for (int i=t;i!=s;i=edge[p[i]].from)
		{
			edge[p[i]].flow+=a[t];
			edge[p[i]^1].flow-=a[t];
		}
		flow+=a[t];
		cost+=(long long)c[t]*a[t];
		return true;
	}
	int solve(int s,int t,int limit_flow,long long &cost)
	{
		int flow=0;
		cost=0;
		while (flow<limit_flow && bellmanford(s,t,limit_flow,flow,cost));
		return flow;
	}
}solver;
struct tree
{
	vector <int> g[maxn];
	int need[maxn],id[maxn],father[maxn];
	int n;
	void init(int n)
	{
		this->n=n;
		memset(need,0,sizeof(need));
		for (int i=0;i<=n;i++)
			g[i].clear();
	}
	void add(int u,int v)
	{
		g[u].push_back(v);
		g[v].push_back(u);
	}
	void update(int k,int x)
	{
		need[k]=x;
	}
	void _dfs(int s)
	{
		dfs(s,-1);
	}
	void dfs(int u,int fa)
	{
		if (fa!=-1) id[u]=id[fa];
		if (need[u]) id[u]=u;
		father[u]=fa;
		for (int i=0;i<g[u].size();i++)
		{
			int v=g[u][i];
			if (v==fa) continue;
			dfs(v,u);
		}
	}
	int _dfs2(int u)
	{
		return dfs2(u,father[u],u);
	}
	int dfs2(int u,int fa,int cur)
	{
		if (u!=cur && need[u]) return need[u];
		int ret=0;
		for (int i=0;i<g[u].size();i++)
		{
			int v=g[u][i];
			if (v==fa) continue;
			ret+=dfs2(v,u,cur);
			
		}
		return ret;
	}
	void debug()
	{
		for (int i=1;i<=n;i++)
			printf("%d ",need[i]);
		printf("\n");
		for (int i=1;i<=n;i++)
			printf("%d ",id[i]);
		printf("\n");
	}
}cand1,cand2;
int a[maxn];
int main()
{
	int n,x,y,q,u,v,k,p;
	scanf("%d%d%d",&n,&x,&y);
	solver.init(2*n+2);
	int s=0,t=2*n+1;
	cand1.init(n);
	cand2.init(n);
	for (int i=1;i<=n;i++)
		scanf("%d",&a[i]);
	for (int i=1;i<n;i++)
	{
		scanf("%d%d",&u,&v);
		cand1.add(u,v);
	}
	for (int i=1;i<n;i++)
	{
		scanf("%d%d",&u,&v);
		cand2.add(u,v);
	}
	scanf("%d",&q);
	for (int i=0;i<q;i++)
	{
		scanf("%d%d",&k,&p);
		cand1.update(k,p);
	}
	scanf("%d",&q);
	for (int i=0;i<q;i++)
	{
		scanf("%d%d",&k,&p);
		cand2.update(k,p);
	}
	cand1._dfs(x);
	cand2._dfs(y);
//	cand1.debug();
//	cand2.debug();
	int flows=0;
	bool flag=true;
	for (int i=1;i<=n;i++)
	{
		if (cand1.need[i])
		{
			solver.addedge(s,i,cand1.need[i]-cand1._dfs2(i),0);
			if (cand1.need[i]-cand1._dfs2(i)<0) flag=false;
			flows+=cand1.need[i]-cand1._dfs2(i);
		}
	}
	int flowt=0;
	for (int i=1;i<=n;i++)
	{
		if (cand2.need[i])
		{
			solver.addedge(i+n,t,cand2.need[i]-cand2._dfs2(i),0);
			if (cand2.need[i]-cand2._dfs2(i)<0) flag=false;
			flowt+=cand2.need[i]-cand2._dfs2(i);
		}
	}
	for (int i=1;i<=n;i++)
		solver.addedge(cand1.id[i],cand2.id[i]+n,1,-a[i]);
	long long cost;
	if (!flag || flows!=flowt || solver.solve(s,t,INF,cost)!=flows)
		printf("-1\n");
	else printf("%I64d\n",-cost);
	return 0;
}
