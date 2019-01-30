#include <bits/stdc++.h>
using namespace std;
typedef unsigned long long ull;
typedef long long ll;
const int hb=61;
ull tr[27],hs[255][255];
ll st[255][255];
char s[255][255];
struct Manacher
{
	ull p[505],s[505];
	int len[505];
	int init(const vector <ull> &p)
	{
		int n=p.size();
		for (int i=1,j=0;i<=2*n;j++,i+=2)
	    {
	        s[i]=-1;
	        s[i+1]=p[j];
	    }
	    s[0]=-2;
	    s[2*n+1]=-1;
	    s[2*n+2]=-3;
	    s[2*n+3]=-4;
	    return 2*n+1;
	}
	void manacher(int n)
	{
	    int mx=0,p=0;
	    for (int i=1;i<=n;i++)
	    {
	        if(mx>i) 
				len[i]=min(mx-i,len[2*p-i]);
	        else len[i]=1;
	        while (s[i-len[i]]==s[i+len[i]])
				len[i]++;
	        if (len[i]+i>mx)
				mx=len[i]+i,p=i;
	    }
	}	
	int solve(const vector<ull> &p)
	{
		if (!p.size()) return 0;
		int n=init(p);
	    memset(len,0,sizeof(len));
	    manacher(n);
	    int ans=0;
	    for (int i=1;i<=n;i++)
			ans+=len[i]/2;
		return ans;
	}
}solver;
ll lowbit(ll x)
{
	return x&(x-1);
}
int main()
{
	tr[0]=1;
	for (int i=1;i<=26;i++)
		tr[i]=tr[i-1]*hb;
	int n,m;
	scanf("%d%d",&n,&m);
	for (int i=0;i<n;i++)
		scanf("%s",&s[i]);
	for (int i=0;i<n;i++)
	{
		hs[i][0]=tr[s[i][0]-'a'+1];
		st[i][0]=(1ll<<(s[i][0]-'a'));
		for (int j=1;j<m;j++)
		{
			hs[i][j]=hs[i][j-1]+tr[s[i][j]-'a'+1];
			st[i][j]=st[i][j-1]^(1ll<<(s[i][j]-'a'));
		}	
	}
	ll ans=0;
	for (int l=0;l<m;l++)
		for (int r=l;r<m;r++)
		{
			vector<ull> p;
			for (int i=0;i<n;i++)
			{
				ull cur=hs[i][r]-(l==0?0:hs[i][l-1]);
				ll tmp=st[i][r]^(l==0?0:st[i][l-1]);
				if (lowbit(tmp)==0)
					p.push_back(cur);
				else
				{
					ans+=solver.solve(p);
					p.clear();
				}				
			}
			ans+=solver.solve(p);
		}
	printf("%lld\n",ans);
	return 0;
}
