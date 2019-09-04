### **ACM图论专题**
#### **1、Tarjan算法**
```C++
struct Edge  
{  
    int u,v;  
    Edge(int u,int v):u(u),v(v){}  
    bool operator <(const Edge &rhs)const  
    {  
        return (u<rhs.u) || (u==rhs.u && v<rhs.v);  
    }  
    bool operator ==(const Edge &rhs)const  
    {  
        return (u==rhs.u) && (v==rhs.v);  
    }  
};  
set <Edge> bridge;  
vector <int> num;  
vector <int> g[maxn];   
int pre[maxn],ebccno[maxn],dfs_cnt,ebcc_cnt;  
bool is_cut[maxn];  
void init(int n)  
{  
    for (int i=0;i<n;i++)  
        g[i].clear();  
    bridge.clear();  
    num.clear();  
}  
int dfs1(int u,int fa)  
{  
    int lowu,child=0;  
    lowu=pre[u]=++dfs_cnt;  
    for (int i=0;i<g[u].size();i++)  
    {  
        int v=g[u][i];  
        if (!pre[v])  
        {  
            child++;  
            int lowv=dfs1(v,u);  
            lowu=min(lowu,lowv);  
            if (lowv>=pre[u])  
            {  
                is_cut[u]=true;  
                if (lowv>pre[u])  
                {  
                    bridge.insert(Edge(u,v));  
                    bridge.insert(Edge(v,u));  
                }     
            }  
        }  
        else if (pre[v]<pre[u] && v!=fa)  
            lowu=min(lowu,pre[v]);  
    }  
    if (fa<0 && child==1) is_cut[u]=false;  
    return lowu;  
}  
void dfs2(int u,int fa)  
{  
    if (bridge.count(Edge(u,fa))) return;  
    ebccno[u]=ebcc_cnt;  
    for (int i=0;i<g[u].size();i++)  
    {  
        int v=g[u][i];  
        if (!ebccno[v]) dfs2(v,u);  
    }  
}  
/*	
	求边双连通分量
*/
void find_ebcc(int n)  
{  
    memset(pre,0,sizeof(pre));  
    memset(is_cut,false,sizeof(is_cut));  
    memset(ebccno,0,sizeof(ebccno));  
    dfs_cnt=ebcc_cnt=0;  
    for (int i=0;i<n;i++)  
        if (!pre[i])  
            dfs1(i,-1);  
    for (int i=0;i<n;i++)  
        if (!ebccno[i])  
        {  
            ebcc_cnt++;  
            dfs2(i,-1);  
        }  
}   
```
#### **2、LCA算法**
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
const int mod=1e9+7;
vector <int> g[100000+5];
const int LN=20;
int par[100000+5][20],lev[100000+5];
int st[100000+5],en[100000+5],dfs_clock=0;
void dfs(int u,int p)
{
	st[u]=++dfs_clock;
    if (p!=-1)
    {
        lev[u]=lev[p]+1;
        par[u][0]=p;
        for (int j=1;j<LN;j++)
            if (par[u][j-1]!=-1)
                par[u][j]=par[par[u][j-1]][j-1];
    }
    for (int i=0;i<g[u].size();i++)
    {
    	int v=g[u][i];
    	if (v==p) continue;
    	dfs(v,u);
	}
	en[u]=++dfs_clock;
}
int lca(int u,int v)
{
    if (lev[u]<lev[v])
        swap(u,v);
    for (int j=LN-1;j>=0;j--)
        if (par[u][j]+1 && lev[par[u][j]]>=lev[v])
            u=par[u][j];
    if (u==v)
        return u;
    for (int j=LN-1;j>=0;j--)
    {
        if (par[u][j]+1 && par[u][j]!=par[v][j])
        {
            u=par[u][j];
            v=par[v][j];
        }
    }
    return par[u][0];
}
int main()
{
	int n,u,v;
	scanf("%d",&n);
	for (int i=1;i<n;i++)
	{
		scanf("%d%d",&u,&v);
		g[u].push_back(v);
		g[v].push_back(u);
	}
	memset(par,-1,sizeof(par));
	dfs(1,-1);
}
```
#### **3、网络流**
##### （1）Dinic
```C++
struct Edge  
{  
    int from,to,cap,flow;  
    Edge (int from,int to,int cap,int flow):from(from),to(to),cap(cap),flow(flow){}  
};  
vector <Edge> edge;  
vector <int> g[maxn];  
int cur[maxn],d[maxn];  
bool vis[maxn];  
void init(int n)  
{  
    edge.clear();  
    for (int i=0;i<=n;i++)  
        g[i].clear();  
}  
void addedge(int from,int to,int cap)  
{  
    edge.push_back(Edge(from,to,cap,0));  
    edge.push_back(Edge(to,from,0,0));  
    int m=edge.size();  
    g[from].push_back(m-2);  
    g[to].push_back(m-1);  
}  
bool bfs(int s,int t)  
{  
    memset(vis,false,sizeof(vis));  
    queue <int> q;  
    q.push(s);  
    vis[s]=true;  
    d[s]=0;  
    while (!q.empty())  
    {  
        int x=q.front();  
        q.pop();  
        for (int i=0;i<g[x].size();i++)  
        {  
            Edge& e=edge[g[x][i]];  
            if (!vis[e.to] && e.cap>e.flow)   
            {  
                vis[e.to]=true;  
                d[e.to]=d[x]+1;  
                q.push(e.to);  
            }  
        }  
    }  
    return vis[t];  
}  
int dfs(int x,int a,int t)  
{  
    if (x==t || a==0) return a;  
    int flow=0,f;  
    for (int& i=cur[x];i<g[x].size();i++)  
    {  
        Edge& e=edge[g[x][i]];  
        if (d[x]+1==d[e.to] && (f=dfs(e.to,min(a,e.cap-e.flow),t))>0)  
        {  
            e.flow+=f;  
            edge[g[x][i]^1].flow-=f;  
            flow+=f;  
            a-=f;  
            if (a==0) break;  
        }  
    }  
    return flow;  
}  
int dinic(int s,int t)  
{  
    int flow=0;  
    while (bfs(s,t))  
    {  
        memset(cur,0,sizeof(cur));  
        flow+=dfs(s,INF,t);  
    }  
    return flow;  
}  
```
##### （2）MCMF
```C++
struct Edge  
{  
    int from,to,cap,flow,cost;  
    Edge (int from,int to,int cap,int flow,int cost):from(from),to(to),cap(cap),flow(flow),cost(cost){}  
};  
vector <Edge> edge;  
vector <int> g[maxn];  
int a[maxn],p[maxn],c[maxn],inq[maxn];  
void init(int n)  
{  
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
bool bellmanford(int s,int t,int limit_flow,int &flow,long long& cost)  
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
            Edge& e=edge[g[x][i]];  
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
int MCMF(int s,int t,int limit_flow,long long& cost)  
{  
    int flow=0;  
    cost=0;  
    while (flow<limit_flow && bellmanford(s,t,limit_flow,flow,cost));  
    return flow;  
}  
```
##### （3）Diijkstra优化的MCMF
```C++
const int N=4005,MAX_V=4005,INF=(int)1e9;
typedef pair<int,int> P;
struct edge{int to,cap,cost,rev;};
int n, k;
int a[N];
struct MincostFlow 
{
    int V; //Please set V!!!!
    vector<edge> G[MAX_V];
    int h[MAX_V];
    int dist[MAX_V];
    int prevv[MAX_V],preve[MAX_V];
    void init(int x) 
	{
        V=x;
        for (int i=0;i<V;i++) G[i].clear();
    }
    void add_edge(int from, int to, int cap, int cost)
	{
        G[from].push_back((edge){to,cap,cost,(int)G[to].size()});
        G[to].push_back((edge){from,0,-cost,(int)G[from].size()-1});
    }
    int min_cost_flow(int s, int t, int f)
	{
        int res = 0;
        fill(h,h+V,0);
        while (f>0){
            priority_queue<P,vector<P>,greater<P> > que;
            fill(dist,dist+V,INF);
            dist[s]=0;
            que.push(P(0,s));
            while (!que.empty())
			{
                P p=que.top();que.pop();
                int v=p.second;
                if (dist[v]<p.first) continue;
                for (int i=0;i<G[v].size();i++)
				{
                    edge &e=G[v][i];
                    if (e.cap>0 && dist[e.to]>dist[v]+e.cost+h[v]-h[e.to])
					{
                        dist[e.to]=dist[v]+e.cost+h[v]-h[e.to];
                        prevv[e.to]=v;
                        preve[e.to]=i;
                        que.push(P(dist[e.to],e.to));
                    }
                }
            }
            if (dist[t]==INF) return -1;
            for (int v=0;v<V;v++) h[v]+=dist[v];
            int d=f;
            for (int v=t;v!=s;v=prevv[v])
			{
                d=min(d,G[prevv[v]][preve[v]].cap);
            }
            f-=d;
            res+=d*h[t];
            for (int v=t;v!=s;v=prevv[v])
			{
                edge &e =G[prevv[v]][preve[v]];
                e.cap-=d;
                G[v][e.rev].cap+=d;
            }
        }
        return res;
    }
} mf;
```
#### **4、无向图欧拉通道（输出方案）**
```C++
#include <bits/stdc++.h>
using namespace std;
struct Edge 
{
	int from,to;
};
vector <int> g[300000+5];
vector <Edge> edges;
bool vis[300000+5];
int deg[300000+5],a[100000+5],b[100000+5];
int cur[300000+5];
void addedge(int from,int to)
{
	edges.push_back(Edge{from,to});
	g[from].push_back(edges.size()-1);
	g[to].push_back(edges.size()-1);
}
vector <int> path;
void dfs(int u)
{
	for (int& i=cur[u];i<g[u].size();i++)
	{
		int j=i;
		if (!vis[g[u][j]])
		{
			vis[g[u][j]]=true;
			int v;
			if (edges[g[u][j]].from==u) v=edges[g[u][j]].to;
			else v=edges[g[u][j]].from;
			dfs(v);
			path.push_back(g[u][j]);
		}
	}
}
int main()
{
	int n;
	scanf("%d",&n);
	for (int i=0;i<n-1;i++)
		scanf("%d",&a[i]);
	for (int i=0;i<n-1;i++)
		scanf("%d",&b[i]);
	vector <int> key;
	bool flag=true;
	for (int i=0;i<n-1;i++)
	{
		key.push_back(a[i]);
		key.push_back(b[i]);
		if (a[i]>b[i]) flag=false;
	}
	sort(key.begin(),key.end());
	int k=unique(key.begin(),key.end())-key.begin();
	for (int i=0;i<n-1;i++)
	{
		a[i]=lower_bound(key.begin(),key.begin()+k,a[i])-key.begin();
		b[i]=lower_bound(key.begin(),key.begin()+k,b[i])-key.begin();
		deg[a[i]]++;
		deg[b[i]]++;
	}
	int cnt=0;
	for (int i=0;i<k;i++)
		if (deg[i]%2) cnt++;
	if (cnt>2) flag=false;  //欧拉通道条件 
	if (flag)
	{
		for (int i=0;i<n-1;i++)
			addedge(a[i],b[i]);   //建边 
		if (cnt==2)
		{
			for (int i=0;i<k;i++)
				if (deg[i]%2) 
				{
					dfs(i); 
					break;
				}
		}
		else 
			dfs(0);
		if (path.size()==n-1)  //判断连通性!! 
		{
			vector <int> ans; //输出方案 
			ans.push_back(edges[path[0]].from);
			ans.push_back(edges[path[0]].to);
			bool flag=true;
			for (int i=1;i<path.size();i++)  
			{
				if (edges[path[i]].from==ans.back())
					ans.push_back(edges[path[i]].to);
				else if (edges[path[i]].to==ans.back())
					ans.push_back(edges[path[i]].from);
				else flag=false;
			}
			if (flag)
			{
				for (int i=0;i<ans.size();i++)
					printf("%d%c",key[ans[i]],i==ans.size()-1?'\n':' ');
			}
			else
			{
				ans.clear();
				ans.push_back(edges[path[0]].to);
				ans.push_back(edges[path[0]].from);
				for (int i=1;i<path.size();i++)  
				{
					if (edges[path[i]].from==ans.back())
						ans.push_back(edges[path[i]].to);
					else ans.push_back(edges[path[i]].from);
				}
				for (int i=0;i<ans.size();i++)
					printf("%d%c",key[ans[i]],i==ans.size()-1?'\n':' ');
			}
		}
		else printf("-1\n");
	} 
	else printf("-1\n");
	return 0;
}
```
#### **5、查分约束系统**
##### （1）建图
 ![1.jpg-63.4kB][1]
##### （2）求解
从1开始使用Bellmanford算法求最短路，若有负环则无解。
输出可行解
增加源点，从源点出发在各顶点的最短路为一组可行解
![2.jpg-42.7kB][2]
   
##### （3）代码

```C++
#include <bits/stdc++.h>
using namespace std;
struct Edge
{
	int u,v,c;
};
vector <Edge> edges;
vector <int> g[100000+5];
void addedge(int from,int to,int c)
{
	edges.push_back(Edge{from,to,c});
	int m=edges.size();
	g[from].push_back(m-1);
}
bool inq[100000+5];
int d[100000+5],cnt[100000+5];
bool bellmanford(int n)
{
	memset(inq,0,sizeof(inq));
	memset(d,0x3f,sizeof(d));
	memset(cnt,0,sizeof(cnt));
	queue <int> q;
	q.push(0);
	d[0]=0;
	inq[0]=true;
	while (!q.empty())
	{
		int x=q.front();
		q.pop();
		inq[x]=false;
		for (int i=0;i<g[x].size();i++)
		{
			Edge&e = edges[g[x][i]];
			if (d[e.v]>d[e.u]+e.c)
			{
				d[e.v]=d[e.u]+e.c;
				if (++cnt[e.v]>n)
					return false;
				inq[e.v]=true;
				q.push(e.v);
			}
		}
	}
	return true;
}
int _time[100000+5];
int main()
{
	int n,cas=0;
	while (cin>>n &&n)
	{
		edges.clear();
		for (int i=0;i<=n;i++)
			g[i].clear();
		for (int i=1;i<=n;i++)
			cin>>_time[i];
		string s;
		int a,b;
		while (cin>>s && s!="#")
		{
			cin>>a>>b;
			if (s=="FAS")
				addedge(a,b,_time[a]);
			else if (s=="FAF")
				addedge(a,b,_time[a]-_time[b]);
			else if (s=="SAF")
				addedge(a,b,-_time[b]);
			else addedge(a,b,0);
		}
		for (int i=1;i<=n;i++)
			addedge(0,i,0);
		printf("Case %d:\n",++cas);
		if (bellmanford(n))
		{
			int l=1,r=1e9;
			while (l<=r)
			{
				int mid=(l+r)>>1;
				for (int i=0;i<edges.size();i++)
					if (edges[i].u==0)
						edges[i].c=mid-_time[edges[i].v];
				bool flag=bellmanford(n);
				for (int i=1;i<=n;i++)
					if (d[i]<0) flag=false;
				if (flag) r=mid-1;
				else l=mid+1;
			}
			for (int i=0;i<edges.size();i++)
				if (edges[i].u==0)
					edges[i].c=l-_time[edges[i].v];
			bellmanford(n);
			for (int i=1;i<=n;i++)
				printf("%d %d\n",i,d[i]);
		}
		else
			printf("impossible\n");
		printf("\n");
	}
	return 0;
}
```
#### **6、支配树（DAG）**
问题：我们有一个有向图(可以有环)，定下了一个节点为起点s。现在我们要求：从起点s出发，走向一个点p的所有路径中，必须要经过的点有哪些{xp}。换言之，删掉{xp}中的任意一个点xpi以及它的入边出边，都会使s无法到达p。
性质：
1、它是一棵树，根节点是我们选定的起点s。
2、对于每个点i，它到根的链上的点集就是对于它的必经点集{xi}。
3、对于每个点i，它是它的支配树上的子树内的点的必经点。
```C++
struct Dominant_Tree
{
	vector <int> g[maxn];  //原图
	vector <int> e[maxn];  //反图
	vector <int> t[maxn];  //支配树
	int deg[maxn];          //原图顶点的度
	int a[maxn];            //topo序 
	int n;                   //顶点数
	int f[maxn][20];       //lca
	int dep[maxn];         //支配树上顶点深度 
	Dominant_Tree()
	{
		memset(deg,0,sizeof(deg));
	}
	Dominant_Tree(int n)
	{
		this->n=n;
		for (int i=0;i<=n;i++)
			g[i].clear(),e[i].clear(),t[i].clear();
		memset(deg,0,sizeof(deg));
	} 
	void init(int n)
	{
		this->n=n;
		for (int i=0;i<=n;i++)
			g[i].clear(),e[i].clear(),t[i].clear();
		memset(deg,0,sizeof(deg));
	}
	void addedge(int u,int v)
	{
		g[u].push_back(v);
		e[v].push_back(u);
		deg[v]++;
	} 
	void work()
	{
		topo();
		build();
	}
	void topo()
	{
		int cnt=0;
		queue<int> q;
		for (int i=0;i<n;i++)
		{
			if (deg[i]==0)
				q.push(i);
		}
		while (!q.empty())
		{
			int u=q.front();
			a[cnt++]=u;
			q.pop();
			for (int i=0;i<g[u].size();i++)
			{
				int v=g[u][i];
				deg[v]--;
				if (deg[v]==0)
					q.push(v);
			}
		}
	}
	int lca(int x,int y)
	{
	    if (dep[x]>dep[y]) swap(x,y);
	    for (int i=19;i>=0;i--) 
			if (dep[y]>dep[x] && dep[f[y][i]]>=dep[x]) 
				y=f[y][i];
	    for (int i=19;i>=0;i--) 
			if (f[x][i]!=f[y][i]) 
				x=f[x][i],y=f[y][i];
	    return x==y?x:f[x][0];
	}
	void build()
	{
		int rt=a[0];
		dep[rt]=1;
		for (int i=1;i<n;i++)
		{
			int u=a[i],fa=-1;
			for (int j=0;j<e[u].size();j++)
			{
				int v=e[u][j];
				fa=(fa==-1?v:lca(fa,v));
			}
			dep[u]=dep[fa]+1;
			f[u][0]=fa;
			t[fa].push_back(u);
			for (int i=1;i<=19;i++) 
				f[u][i]=f[f[u][i-1]][i-1];
		}
	}
	void dubug()
	{
		for (int i=0;i<n;i++)
			cout<<f[i][0]<<' ';
		cout<<endl;
	}
	int query(int u,int v)
	{
		return dep[u]+dep[v]-dep[lca(u,v)];
	}
}solver;
```
#### **7、无向图最大团**
给定无向图G=(V,E)。如果U包含于V，且对任意u，v属于U且有(u，v)属于E，则称U是G的完全子图。
G的完全子图U是G的团当且仅当U不包含在G的更大的完全子图中，即U就是最大完全子图。
G的最大团是指G中所含顶点数最多的团。
```C++
#include <bits/stdc++.h>
using namespace std;
const int maxn=1010;
int best;
int num[maxn];
int path[maxn];
int g[maxn][maxn],n;
bool dfs(int *adj,int total,int cnt)
{
    int i,j,k;
    int t[maxn];
    if (total==0)
	{ 
        if (best<cnt)
		{
            best=cnt; 
			return true;
        }
        return false;
    }
    for (i=0;i<total;i++)
	{
		//剪枝1,若当前顶点数量cnt加上还能够增加的最大数量仍小于 best则退出并返回false
        if (cnt+(total-i)<=best) return false;
		// 剪枝2,若当前顶点数量cnt加上包含adj[i]的最大团顶点数仍小于 best则退出并返回false
        if (cnt+num[adj[i]]<=best) return false;
		// 扫描与u相连的顶点中与 adj[u]相连的顶点并存储到数组 t[]中，数量为k 
        for (k=0,j=i+1;j<total;j++) 
            if (g[adj[i]][adj[j]])
                t[k++]=adj[j];
                if (dfs(t,k,cnt+1)) return true;
    } 
	return false;
}
int MaximumClique()
{
    int i, j, k;
    int adj[maxn];
    if (n<=0) return 0;
    best=0;
    for (i=n-1;i>=0;i--)
	{
        for (k=0,j=i+1;j<n;j++)
           if (g[i][j]) adj[k++]=j;
        dfs(adj,k,1);
        num[i]=best;
    }
    return best;
}
int main()
{
    cin>>n;
    for(int i=0;i<n;i++)
    	for(int j=0;j<n;j++)
    		scanf("%d",&g[i][j]);
    int t=MaximumClique();
    cout<<t<<endl;
    return 0;
}
```
#### **8、无向图第k小团**
定义：给定无向图G=(V,E)。如果U包含于V，且对任意u，v属于U且有(u，v)属于E，则称U是G的完全子图。
第k小团定义：给定无向图G=(V,E)，且每个点有对应的权值，无向图的一个团U的权值定义为U中所有点权值之和。第k小团为无向图G中所有团中权值第k小团的权值。
算法：搜索+暴力
分析：
若当前团为U，则点v可以加入U当且仅当v到U中所有点都存在边。
记adj[i]表示i能到达的邻接点的集合，state表示U中点的集合，
则v能加入U当且仅当(adj[v]&state)==state。
计算第k小可以二分答案或者优先级队列。
```C++
#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int n,k;
__int128_t adj[105];
__int128_t bit[105];
ll v[105];
string s[105];
struct node
{
	ll val;
	__int128_t state;
	int pre;
	bool operator < (const node &rhs) const
	{
		return val>rhs.val;
	}
};
int main()
{
	bit[0]=1;
	for (int i=1;i<105;i++)
		bit[i]=bit[i-1]*2;
	cin>>n>>k;
	for (int i=0;i<n;i++)
		cin>>v[i];
	for (int i=0;i<n;i++)
		cin>>s[i];
	for (int i=0;i<n;i++)
		for (int j=0;j<n;j++)
		{
			if (i==j) continue;
			if (s[i][j]=='1') adj[i]+=bit[j];
		} 
	__int128_t tmp=0;
	priority_queue <node> pq;
	pq.push(node{0,tmp,0});
	int cnt=0;
	while (!pq.empty())
	{
		node cur=pq.top();
		pq.pop();
		cnt++;
		if (cnt==k)
		{
			cout<<cur.val<<endl;
			return 0;
		}
		__int128_t state=cur.state;
		int pre=cur.pre;
		for (int i=pre;i<n;i++)
		{
			if ((adj[i]&state)==state)
			{
				state+=bit[i];
				pq.push(node{cur.val+v[i],state,i+1});
				state-=bit[i];
			}
		}
	}
	cout<<-1<<endl;
	return 0;
}
```
#### **9、强联通+缩点**
给定一个图，添加最少的边，使之成为一个强连通图。
分析：先求出强连通分量，缩点之后求出入度和出度为0的点的个数，入度为0的个
数和出度为0的个数的最大值即为答案。
```C++
#include <bits/stdc++.h>
using namespace std;
vector <int> g[20000+5];
vector <int> rg[20000+5];
bool vis[20000+5];
vector <int> p;
int u[50000+5],v[50000+5],in[20000+5],out[20000+5],bcc[20000+5];
void dfs(int u)
{
	if (vis[u]) return;
	vis[u]=true;
	for (int i=0;i<g[u].size();i++)
	{
		int v=g[u][i];
		dfs(v);
	}
	p.push_back(u);
} 
void rdfs(int cnt,int u)
{
	if (vis[u]) return;
	bcc[u]=cnt;
	vis[u]=true;
	for (int i=0;i<rg[u].size();i++)
	{
		int v=rg[u][i];
		rdfs(cnt,v);
	}
}
int main()
{
	int n,m,t;
	cin>>t;
	while (t--)
	{
		cin>>n>>m; 
		for (int i=1;i<=n;i++)
		{
			g[i].clear();
			rg[i].clear();
		}
		p.clear();
		for (int i=0;i<m;i++)
		{
			cin>>u[i]>>v[i];
			g[u[i]].push_back(v[i]);
			rg[v[i]].push_back(u[i]);
		}
		memset(vis,false,sizeof(vis));
		for (int i=1;i<=n;i++)
		{
			if (!vis[i])
				dfs(i);
		}
		reverse(p.begin(),p.end());
		memset(vis,false,sizeof(vis));
		int cnt=0;
		for (int i=0;i<p.size();i++)
		{
			if (!vis[p[i]])
			{
				cnt++;
				rdfs(cnt,p[i]);
			} 
		}
		if (cnt==1)
			cout<<0<<endl;
		else
		{
			for (int i=1;i<=cnt;i++)
				in[i]=out[i]=0;
			for (int i=0;i<m;i++)
			{
				if (bcc[u[i]]!=bcc[v[i]])
					in[bcc[v[i]]]++,out[bcc[u[i]]]++;
			}
			int a=0,b=0;
			for (int i=1;i<=cnt;i++)
			{
				if (in[i]==0) a++;
				if (out[i]==0) b++;
			}
			cout<<max(a,b)<<endl;
		}
	}
	return 0;
}
```


  [1]: http://static.zybuluo.com/chenyuqi/nx7rki8krntinka9a2tsqlnm/1.jpg
  [2]: http://static.zybuluo.com/chenyuqi/p0hep63aheifuzo08nao7g0s/2.jpg
