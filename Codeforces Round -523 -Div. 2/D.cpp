#include <bits/stdc++.h>
using namespace std;
const int maxn=1e5+5;
const int mod=1e9+7;
const int INF=2e9;
struct node
{
	int st,en;
	node(int st=0,int en=0):st(st),en(en){}
	bool operator < (const node &rhs) const
	{
		return (st<rhs.st) || (st==rhs.st && en<rhs.en); 
	}
};
node a[maxn];
int main()
{
	int n,x,y;
	scanf("%d%d%d",&n,&x,&y);
	int k=x/y;
	bool flag=(x%y==0);
	multiset <int> s;
	long long ans=0;
	s.insert(0);
	s.insert(INF);
	for (int i=0;i<n;i++)
		scanf("%d%d",&a[i].st,&a[i].en);
	sort(a,a+n);
	for (int i=0;i<n;i++)
	{
		set<int>::iterator it=s.lower_bound(a[i].st);
		int p=*(--it);
		if (p==0)
		{
			ans=(ans+1ll*x+1ll*y*1ll*(a[i].en-a[i].st))%mod;
			s.insert(a[i].en);
		}
		else
		{
			if ((a[i].st-p<k) || ((a[i].st-p)==k && !flag))
			{
				ans=(ans+1ll*y*1ll*(a[i].en-p))%mod;
				s.erase(s.lower_bound(p));
				s.insert(a[i].en);
			}
			else
			{
				ans=(ans+1ll*x+1ll*y*1ll*(a[i].en-a[i].st))%mod;
				s.insert(a[i].en);
			}
		}
	}
	printf("%d\n",(int)ans);
	return 0;
}
