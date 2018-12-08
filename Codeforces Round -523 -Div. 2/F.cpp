#include <bits/stdc++.h>
using namespace std;
const int maxn=1505;
struct node
{
	int a,b,c;
	node(int a=0,int b=0,int c=0):a(a),b(b),c(c){}
	bool operator < (const node &rhs) const
	{
		return a<rhs.a || (a==rhs.a && b<rhs.b) || (a==rhs.a && b==rhs.b && c<rhs.c);
	}
	bool operator == (const node &rhs) const
	{
		return a==rhs.a && b==rhs.b && c==rhs.c;
	}
};
map<node,bool> query;
int leave[maxn];
//leave[i]==1:is a leave
//leave[i]==0:not checked
//leave[i]==-1:not a leave
bool expect[maxn];
//expect[i]==false:not expected in the random()
//expect[i]==true:expected in the random()
int n,len;
bool interact(int a,int b,int c)
{
	if (query.count(node(a,b,c))) return query[node(a,b,c)];
	if (query.count(node(c,b,a))) return query[node(c,b,a)];
	string respond;
	cout<<"? "<<a<<' '<<b<<' '<<c<<endl;
	fflush(stdout);
	cin>>respond;
	bool flag=(respond=="Yes")?true:false;
	query[node(a,b,c)]=flag;
	return flag;
}
int random()
{
	int cnt=0;
	while (++cnt<n)
	{
		int k=rand()%n+1;
		if (expect[k]) return k; 
	}
	for (int i=1;i<=n;i++)
		if (expect[i]) return i;
	assert(false);
}
bool check_leave(int k)
{
//	cout<<"check_leave: "<<k<<endl;
	memset(expect,true,sizeof(expect));
	expect[k]=false;
	int p=random();
	expect[p]=false;
	int cnt=0;
	while ((cnt++)<n-2)
	{
		int q=random();
		expect[q]=false;
		if (interact(p,k,q)) return false;
	}
//	cout<<"check_leave true"<<endl;
	return true;
}
int get_leave()
{
	while (true)
	{
		int k=rand()%n+1;
//		cout<<"get_leave: "<<k<<endl;
		if (leave[k]==0 && check_leave(k)) 
		{
			leave[k]=1;
			return k; 
		}
		leave[k]=-1;
	}
	assert(false);
}
bool check_leaves(int p,int q)
{
	assert(p!=q);
	int cnt=0;
	for (int i=1;i<=n;i++)
	{
		if (i==p || i==q) continue;
		if (interact(p,i,q))
			cnt++;
	}
//	cout<<"check_leaves: "<<p<<' '<<q<<' '<<cnt<<endl;
	return cnt==len;
}
int get_ano_leave(int p)
{
//	cout<<"get_ano_leave: "<<p<<endl;
	memset(expect,true,sizeof(true));
	for (int i=1;i<=n;i++)
		if (leave[i]!=0)
			expect[i]=false;
	while (true)
	{
		int k=rand()%n+1;
		if (leave[k]==0)
		{
			if (check_leave(k))
			{
				leave[k]=1;
				if (check_leaves(p,k))
					return k;	
			}
			else leave[k]=-1;
		}
	}
}
vector <int> path;
int get_path(int p,int q)
{
//	cout<<"get_path: "<<p<<' '<<q<<endl;
	assert(p!=q);
	for (int i=1;i<=n;i++)
	{
		if (i==p || i==q) continue;
		if (interact(p,i,q))
			path.push_back(i);
	}
} 
int check_root(int r,int p,int q)
{
//	cout<<"check_root: "<<r<<' '<<p<<' '<<q<<endl;
	assert(p!=q);
	assert(p!=r);
	assert(q!=r);
	int cntp=0,cntq=0;
	for (int i=1;i<=n;i++)
	{
		if (i==p || i==r) continue;
		if (interact(p,r,i)) cntp++;
	}
	for (int i=1;i<=n;i++)
	{
		if (i==q || i==r) continue;
		if (interact(q,r,i)) cntq++;
	}
	return cntp==cntq;
}
int get_root(int p,int q)
{
//	cout<<"get root: "<<p<<' '<<q<<endl;
	assert(path.size());
	memset(expect,false,sizeof(expect));
	for (int i=0;i<path.size();i++)
		expect[path[i]]=true;
	while (true)
	{
		int k=random();
		if (check_root(k,p,q)) 
			return k;
	}
	assert(false);
}
void solve()
{
	int p=get_leave();
	int q=get_ano_leave(p);
	get_path(p,q);
	int r=get_root(p,q);
	cout<<"! "<<r<<endl;	
}
int main()
{
	int k;
	cin>>n>>k;
	len=0;
	int p=1,s=0;
	while (true)
	{
		s+=p;
		p*=k;
		len++;
		if (s==n) break;
	}
	len=2*len-3;
//	cout<<len<<endl;
	solve();
	return 0;
}
