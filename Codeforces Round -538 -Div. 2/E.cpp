#include <bits/stdc++.h>
using namespace std;
int gcd(int a,int b)
{
	return b==0?a:gcd(b,a%b);
}
int Rand()
{
	return (rand()<<15) | rand();
	//rand()×î´óÖµ32767 
}
bool vis[1000000+5];
vector <int> v;
int main()
{
	srand(6666);
	int n,res;
	cin>>n;
	int l=n-1,r=1e9;
	int cnt=0;
	while (l<=r)
	{
		cnt++;
		int mid=(l+r)/2;
		cout<<"> "<<mid<<endl;
		fflush(stdout);
		cin>>res;
		if (res==1)
			l=mid+1;
		else r=mid-1; 
	}
	int maxn=l;
	v.push_back(l);
	memset(vis,false,sizeof(vis));
	while (cnt++<60)
	{
		bool flag=true;
		for (int i=1;i<=n && flag;i++)
			if (!vis[i])
				flag=false;
		if (flag) break;
		int k;
		while (true)
		{
			k=Rand()%n+1;
			if (!vis[k])  break;
		}
		vis[k]=true;
		cout<<"? "<<k<<endl;
		fflush(stdout);
		cin>>res;
		v.push_back(res);
	}
	sort(v.begin(),v.end());
	int sz=unique(v.begin(),v.end())-v.begin();
	int d=-1;
	for (int i=0;i<sz;i++)
		for (int j=i+1;j<sz;j++)
		{
			if (d==-1) d=v[j]-v[i];
			else d=gcd(d,v[j]-v[i]);
		}
	int x1=maxn-(n-1)*d;
	cout<<"! "<<x1<<' '<<d<<endl;
	fflush(stdout);
	return 0;
}
