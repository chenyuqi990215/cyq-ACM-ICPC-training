#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
int m,n,f[1<<21];
ll bit[55];
char s[55];
int main()
{
	int m;
	scanf("%d%d",&m,&n);
	vector <int> v;
	map <string,int> id;
	int k=0,op;
	for (int i=0;i<n;i++)
		bit[i]=1ll<<i;
	for (int i=0;i<=m;i++)
	{
		if (i<m) scanf("%d",&op);
		if (op==1 || i==m)
		{
			int sz=unique(v.begin(),v.end())-v.begin();
			ll tmp=0;
			for (int i=0;i<sz;i++)
				tmp|=1ll<<v[i];
			for (int i=0;i<sz;i++)
				bit[v[i]]|=tmp;
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
	int n1=n/2,n2=n-n1;
	for (int i=0;i<n1;i++)
		f[1<<i]=1;
	for (int i=0;i<(1<<n1);i++)
		for (int j=0;j<n1;j++)
			if (i&(1<<j))
				f[i]=max(f[i],f[i&(~bit[j])]+1);
	int ans=0;
	for (int i=0;i<(1<<n2);i++)
	{
		bool flag=true;
		for (int j=0;j<n2;j++)
			if (i&(1<<j)) 
				if ((bit[j+n1]&((ll)i<<n1))!=(1ll<<(j+n1)))
					flag=false;
		if (flag)
		{
			ll tmp=0;
			int cnt=0;
			for (int j=0;j<n2;j++)
				if (i&(1<<j))
					tmp|=bit[j+n1],cnt++;
			int st=(tmp&((1ll<<n1)-1))^((1ll<<n1)-1);
			ans=max(ans,cnt+f[st]);
		}
	}
	printf("%d\n",ans);
	return 0;
} 
