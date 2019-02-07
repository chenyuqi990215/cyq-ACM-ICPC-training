/*
	无向图的最大团==该无向图补图的最大独立集
*/
#include <bits/stdc++.h>
using namespace std;
bool a[55][55];//a为图的邻接表(从1开始) 
int ans,cnt[55],group[55],n,vis[55];
bool dfs(int u,int pos) 
{
    int i,j;
    for (i=u+1;i<=n;i++)
	{
        if (cnt[i]+pos<=ans) return 0; 
        if (a[u][i]) 
		{ 
            for (j=0;j<pos;j++) 
				if (!a[i][vis[j]]) break; 
            if (j==pos)
			{   
                vis[pos]=i;
                if (dfs(i,pos+1)) return 1;    
            }    
        }
    }    
    if (pos>ans)
	{
        for( i=0;i<pos;i++)
            group[i]=vis[i];
        ans=pos;
        return 1;    
    }    
    return 0;
} 
int maxclique()
{
    ans=-1;
    for (int i=n;i>0;i--)
    {
        vis[0]=i;
        dfs(i,1);
        cnt[i]=ans;
    }
    if (ans<0) ans=0;
    return ans;
}
char s[55];
int main()
{
	int m;
	scanf("%d%d",&m,&n);
	memset(a,0,sizeof(a));
	vector <int> v;
	map <string,int> id;
	int k=1,op;
	for (int i=0;i<=m;i++)
	{
		if (i<m) scanf("%d",&op);
		if (op==1 || i==m)
		{
			int sz=unique(v.begin(),v.end())-v.begin();
			for (int i=0;i<sz;i++)
				for (int j=i+1;j<sz;j++)
					a[v[i]][v[j]]=a[v[j]][v[i]]=1;
			v.clear();
		}
		else
		{
			scanf("%s",s);
			if (!id.count(s))
				id[s]=k++;
			v.push_back(id[s]);
		}
	}
	for (int i=1;i<=n;i++)
		for (int j=1;j<=n;j++)
			a[i][j]=!a[i][j];
	printf("%d\n",maxclique());
	return 0;
}
