#include <bits/stdc++.h>
using namespace std;
const int maxb=20;
int a[500000+5],b[maxb+5],p[maxb+5];
void insert(int k)
{
	int x=a[k];
	for (int j=maxb;j>=0;j--)
	{
		if ((x>>j)&1)
		{
			if (b[j]==0)   //如果这一维度是0，则保存当前数作为基中的数 
			{
				b[j]=x;
				p[j]=k;   //保存位置 
				break;
			}
			if (p[j]<k)  //保存最后面的解 
			{
				swap(b[j],x);
				swap(p[j],k);
				//贪心，将x插入这个位置，然后尝试将原来的数继续插入线性基中 
			}
			x^=b[j];
		}
	}
}
int querymax(int l)
{
	int ans=0;
	for (int i=maxb;i>=0;i--)
		if (p[i]>=l && ((ans^b[i])>ans))
			ans^=b[i];
	return ans;
}
struct node
{
	int l,r,id;
	bool operator < (const node &rhs) const
	{
		return r<rhs.r;
	}
}q[500000+5];
int ans[500000+5];
int main()
{
	int n,m;
	scanf("%d",&n);
	for (int i=1;i<=n;i++)
		scanf("%d",&a[i]);
	scanf("%d",&m);
	for (int i=0;i<m;i++)
	{
		scanf("%d%d",&q[i].l,&q[i].r);
		q[i].id=i;
	}
	sort(q,q+m);
	int r=1;
	for (int i=0;i<m;i++)
	{
		while (r<=q[i].r)
			insert(r++);
		ans[q[i].id]=querymax(q[i].l);
	}
	for (int i=0;i<m;i++)
		printf("%d\n",ans[i]);
	return 0;
}
