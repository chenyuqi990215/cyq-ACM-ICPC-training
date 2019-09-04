### **ACM数据结构专题**
#### **1、带权并查集**
```C++
int father(int x)
{
	if (fa[x]==x)
		return x;
	else
	{
		int y=father(fa[x]);
		d[x]+=d[fa[x]];
		fa[x]=y;
		return fa[x];
	}
}
void unionn(int from,int to)
{
	int r1=father(from);
	int r2=father(to);
	if (r1!=r2)
	{
		fa[r1]=r2;
		d[r1]=cnt[r2];
		cnt[r2]+=cnt[r1];
	}
}
```
#### **2、线段树**
##### （1）区间修改 区间累乘
```C++
#include <bits/stdc++.h>
#define mp make_pair
#define lc 2*o
#define rc 2*o+1
#define mc l+(r-l)/2
using namespace std;
typedef long long ll;
typedef pair<int,int> pii;
struct segmul
{
	ll sumv[maxn],mulv[maxn],v;
	int n,Y1,Y2;
	void pushup(int o)
	{
		sumv[o]=sumv[lc]*sumv[rc]%mod;
	}
	void pushdown(int o,int m)
	{
		mulv[lc]=(mulv[lc]*mulv[o])%mod;
		mulv[rc]=(mulv[rc]*mulv[o])%mod;
		sumv[lc]=(sumv[lc]*pow_mod(mulv[o],(m-(m>>1))))%mod;
		sumv[rc]=(sumv[rc]*pow_mod(mulv[o],(m>>1)))%mod;
		mulv[o]=1;
	}
	void _update(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
		{
			sumv[o]=(sumv[o]*pow_mod(v,(r-l+1)))%mod;
			mulv[o]=(mulv[o]*v)%mod;
		}
		else
		{
			pushdown(o,r-l+1);
			int mc=(l+r)>>1;
			if (Y1<=mc) _update(lc,l,mc);
			if (Y2>mc) _update(rc,mc+1,r);
			pushup(o);
		}	
	}
	ll _query(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
			return sumv[o];
		else
		{
			pushdown(o,r-l+1);
			int mc=(l+r)>>1;
			ll ans=1;
			if (Y1<=mc) ans=(ans*_query(lc,l,mc))%mod;
			if (Y2>mc) ans=(ans*_query(rc,mc+1,r))%mod;
			return ans;
		}
	}
	void update(int Y1,int Y2,ll v)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		this->v=v;
		_update(1,1,n); 
	}
	ll query(int Y1,int Y2)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		return _query(1,1,n);
	}
}trmul;
```
##### （2）区间修改 区间求or和
```C++
struct segbit
{
	bitset<65> bitv[maxn],sbitv[maxn],v;
	int n,Y1,Y2;
	void pushup(int o)
	{
		bitv[o]=bitv[lc]|bitv[rc];
	}
	void pushdown(int o)
	{
		bitv[lc]|=sbitv[o];
		bitv[rc]|=sbitv[o];
		sbitv[lc]|=sbitv[o];
		sbitv[rc]|=sbitv[o];
		sbitv[o]=0;
	}
	void _update(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
		{
			bitv[o]|=v;
			sbitv[o]|=v;
		}
		else
		{
			pushdown(o);
			int mc=(l+r)>>1;
			if (Y1<=mc) _update(lc,l,mc);
			if (Y2>mc) _update(rc,mc+1,r);
			pushup(o);
		}	
	}
	bitset<65> _query(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
			return bitv[o];
		else
		{
			pushdown(o);
			int mc=(l+r)>>1;
			bitset<65> ans=0;
			if (Y1<=mc) ans|=_query(lc,l,mc);
			if (Y2>mc) ans|=_query(rc,mc+1,r);
			return ans;
		}
	}
	void update(int Y1,int Y2,bitset<65> v)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		this->v=v;
		_update(1,1,n); 
	}
	bitset<65> query(int Y1,int Y2)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		return _query(1,1,n);
	}
}trbit;
```
##### （3）单点修改，区间查询
```C++
struct Segtree
{
	pii maxv[4*maxn],v;
	int p,ql,qr,n;
	void init(int n)
	{
		this->n=n;
		for (int i=0;i<=4*n;i++)
			maxv[i]=mp(-1,-1);
	}
	void update(int o,int l,int r)
	{
		if (l==r)
			maxv[o]=max(maxv[o],v);
		else
		{
			if (p<=mc) update(lc,l,mc);
			else update(rc,mc+1,r);
			maxv[o]=max(maxv[lc],maxv[rc]);
		}
	}
	pii quiry(int o,int l,int r)
	{
		pii ans=mp(-1,-1);
		if (ql<=l && r<=qr) return maxv[o];
		if (ql<=mc) ans=max(ans,quiry(lc,l,mc));
		if (qr>mc) ans=max(ans,quiry(rc,mc+1,r));
		return ans;
	}
	void update(int p,pii v)
	{
		this->p=p;
		this->v=v;
		update(1,1,n);
	}
	pii quiry(int ql,int qr)
	{
		this->ql=ql;
		this->qr=qr;
		return quiry(1,1,n);
	}
}tr;
```
##### （4）区间加，查询max
```C++
#include <bits/stdc++.h>
#define mp make_pair
#define lc 2*o
#define rc 2*o+1
using namespace std;
typedef long long ll;
const int maxn=4e5+20;
struct segadd
{
	ll maxv[maxn],addv[maxn],v;
	int n,Y1,Y2;
	void init(int n)
	{
		this->n=n;
		for (int i=0;i<4*n;i++)
			maxv[i]=addv[i]=0;
	}
	void pushup(int o)
	{
		maxv[o]=max(maxv[lc],maxv[rc]);
	}
	void pushdown(int o,int m)
	{
		addv[lc]+=addv[o];
		addv[rc]+=addv[o];
		maxv[lc]+=addv[o];
		maxv[rc]+=addv[o];
		addv[o]=0;
	}
	void _update(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
		{
			maxv[o]+=v;
			addv[o]+=v;
		}
		else
		{
			pushdown(o,r-l+1);
			int mc=(l+r)>>1;
			if (Y1<=mc) _update(lc,l,mc);
			if (Y2>mc) _update(rc,mc+1,r);
			pushup(o);
		}	
	}
	ll _query(int o,int l,int r)
	{
		if (Y1<=l && r<=Y2)
			return maxv[o];
		else
		{
			pushdown(o,r-l+1);
			int mc=(l+r)>>1;
			ll ans=-1e16;   //!!!
			if (Y1<=mc) ans=max(ans,_query(lc,l,mc));
			if (Y2>mc) ans=max(ans,_query(rc,mc+1,r));
			return ans;
		}
	}
	void update(int Y1,int Y2,ll v)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		this->v=v;
		_update(1,1,n); 
	}
	ll query(int Y1,int Y2)
	{
		this->Y1=Y1;
		this->Y2=Y2;
		return _query(1,1,n);
	}
}tradd;
```
##### （5） 题意：给出一个长度为n初值为0的数组，以及长度为n的b数组，然后q次操作，add(l,r) 使得区间l~r所有元素+1，或者查询l~r区间$\frac{a[i]}{b[i]} $的和。
分析：考虑到$\frac{a[i]}{b[i]} $最多更新nlogn次，所以每次需要更新时，暴力到叶子节点，对于每一个没叶子节点，维护一个懒标记，表示当前节点至少还要操作几次，才会有节点需要被更新，注意到$ tagv[o]=min(tagv[lc],tagv[rc); $。
```C++
#include <iostream>  
#include <cstring>  
#define maxn 100000+5  
using namespace std;  
int sumv[4*maxn],addv[4*maxn],tagv[4*maxn],b[maxn];  
void init()  
{  
    memset(sumv,0,sizeof(sumv));  
    memset(addv,0,sizeof(addv));  
    memset(tagv,0,sizeof(tagv));  
}  
void push_down(int o)  
{  
    int lc=2*o,rc=2*o+1;  
    if (addv[o])  
    {  
        int v=addv[o];  
        addv[o]=0;  
        addv[lc]+=v;  
        addv[rc]+=v;  
        tagv[lc]-=v;  
        tagv[rc]-=v;  
    }  
}  
void push_up(int o)  
{  
    int lc=2*o,rc=2*o+1;  
    sumv[o]=sumv[lc]+sumv[rc];  
    tagv[o]=min(tagv[lc],tagv[rc]);  
}  
int y1,y2;  
void update(int o,int l,int r)  
{  
    int lc=2*o,rc=2*o+1;  
    if (y1<=l && r<=y2)  
    {  
        tagv[o]--;  
        if (tagv[o]>0)  
        {  
            addv[o]++;  
            return;  
        }  
        if (l==r && tagv[o]<=0)  
        {  
            tagv[o]=b[l];  
            sumv[o]++;  
            return;  
        }  
    }  
    push_down(o);  
    int m=l+(r-l)/2;  
    if (y1<=m) update(lc,l,m);  
    if (y2>m) update(rc,m+1,r);  
    push_up(o);  
}  
int _sum=0;  
void query(int o,int l,int r)  
{  
    int lc=2*o,rc=2*o+1;  
    if (y1<=l && r<=y2)  
    {  
        _sum+=sumv[o];  
        return;  
    }  
    int m=l+(r-l)/2;  
    if (y1<=m) query(lc,l,m);  
    if (y2>m) query(rc,m+1,r);  
}  
void build(int o,int l,int r)  
{  
    int lc=2*o,rc=2*o+1;  
    if (l==r)  
    {  
        tagv[o]=b[l];  
        return;  
    }  
    int m=l+(r-l)/2;  
    build(lc,l,m);  
    build(rc,m+1,r);  
    push_up(o);  
}  
int main()  
{  
    int n,q;  
    char s[10];  
    while (scanf("%d%d",&n,&q)==2)  
    {  
        init();  
        for (int i=1;i<=n;i++)  
            scanf("%d",&b[i]);  
        build(1,1,n);  
        for (int i=0;i<q;i++)  
        {  
            scanf("%s%d%d",s,&y1,&y2);  
            if (strcmp(s,"add")==0)  
                update(1,1,n);  
            else  
            {  
                _sum=0;  
                query(1,1,n);  
                printf("%d\n",_sum);  
            }  
        }  
    }  
    return 0;  
}  
```
3、主席树
##### （1）动态第k大
```C++
#include <cstdio>
#include <cstring>
#include <vector>
#include <algorithm>
using namespace std;
const int maxn=1e5+100;
int a[maxn];
int rt[maxn];//rt[i]表示由数组前i个元素组成的线段树的根结点
struct node
{
    int l,r;//线段树左右子结点点
    int sum;//结点信息，表示这颗子树存在的元素的数目
}T[maxn*20];
int tot=0;//结点编号
vector<int> v;
int getid(int k)
{
    return lower_bound(v.begin(),v.end(),k)-v.begin()+1;
}
void build(int &o,int l,int r)//建立一颗空树
{
    o=++tot;
    T[o].sum=0;
    if(l==r) return;
    int mid=(l+r)/2;
    build(T[o].l,l,mid);
    build(T[o].r,mid+1,r);
}
void update(int l,int r,int &now,int last,int k)
{
    T[++tot]=T[last];//复制线段树
    //更新当前线段树的根结点
    now=tot;
    T[tot].sum++;
    if(l==r) return;//修改到叶子结点为止
    //根据需要修改的k来确定是修改左子树还是修改右子树
    int mid=(l+r)/2;
    if(k<=mid)
        update(l,mid,T[now].l,T[last].l,k);
    else
        update(mid+1,r,T[now].r,T[last].r,k);
}
int query(int l,int r,int x,int y,int k)//查询区间【x，y】中第小的数
{
    if(l==r) return l;//查询到叶子结点为止
    int mid=(l+r)/2;
    int cnt=T[T[y].l].sum-T[T[x].l].sum;//第y颗树比第x颗树在左子树上多的结点数
    if(cnt>=k)//答案在左子树上
        return query(l,mid,T[x].l,T[y].l,k);
    else
        return query(mid+1,r,T[x].r,T[y].r,k-cnt);
}
int main()
{
    int n,m;
    while(~scanf("%d%d",&n,&m))
    {
        for(int i=1;i<=n;i++)
        {
            scanf("%d",&a[i]);
            v.push_back(a[i]);
        }
        sort(v.begin(),v.end());
        v.erase(unique(v.begin(),v.end()),v.end());
        build(rt[0],1,n);
        for(int i=1;i<=n;i++)
            update(1,n,rt[i],rt[i-1],getid(a[i]));
        while(m--)
        {
            int x,y,k;
            scanf("%d%d%d",&x,&y,&k);
            printf("%d\n",v[query(1,n,rt[x-1],rt[y],k)-1]);
        }
    }
    return 0;
}
```
##### （2）主席树查询区间小于等于h的个数
```C++
#include <bits/stdc++.h>
using namespace std;
const int maxn=1e5+100;
vector<int> v;
int n,m,a[maxn],rt[maxn],tot;
struct node
{
    int l,r,sum;
}T[maxn*20];
void build(int &o,int l,int r)
{
    o=++tot;
    T[o].sum=0;
    if (l==r) return;
    int mid=(l+r)>>1;
    build(T[o].l,l,mid);
    build(T[o].r,mid+1,r);
}
void update(int l,int r,int &now,int last,int k)
{
    T[++tot]=T[last];
    now=tot;
    T[now].sum++;
    if (l==r) return;
    int mid=(l+r)>>1;
    if (k<=mid) update(l,mid,T[now].l,T[last].l,k);
    else update(mid+1,r,T[now].r,T[last].r,k);
}
int query(int l,int r,int x,int y,int k)
{
    if (l==r) return T[y].sum-T[x].sum;
    int mid=(l+r)>>1;
    if (k<=mid) return query(l,mid,T[x].l,T[y].l,k);
    else
    {
        int ret=0;
        ret+=T[T[y].l].sum-T[T[x].l].sum;
        ret+=query(mid+1,r,T[x].r,T[y].r,k);
        return ret;
    }
}
int main()
{
	int t;
    scanf("%d",&t);
    for(int kase=1;kase<=t;kase++)
{
tot=0;   //!!!
        scanf("%d%d",&n,&m);
        for (int i=1;i<=n;i++) scanf("%d",&a[i]);
        build(rt[0],1,n);
        for(int i=1;i<=n;i++)
            update(1,n,rt[i],rt[i-1],a[i]);
        for(int i=1;i<=m;i++)
        {
            int l,r,h;
            scanf("%d%d%d",&l,&r,&h);
            printf("%d\n",query(1,n,rt[l-1],rt[r],h));
        }
    }
    return 0;
}
```
#### **4、RMQ**
```C++
void get_RMQ(int r)  
{  
    for (int i=0;i<r;i++)  
        f[0][i]=times[i];  
    for (int i=1;(1<<i)<r;i++)  
        for (int j=0;(j+(1<<i))<=r;j++)  
            f[i][j]=max(f[i-1][j],f[i-1][j+(1<<(i-1))]);    
}  
int quiry_RMQ(int l,int r)  
{  
    int k=0;  
    while ((1<<(k+1))<=(r-l+1)) k++;  
    return max(f[k][l],f[k][r-(1<<k)+1]);  
}  
```
#### **5、分块算法**
**基本代码**
```C++
//num:分块的个数  
//belong[i]:i所属的块  
//block:块的大小  
//l[i]:i这块的左端点位置  
//r[i]:i这块的右端点位置   
//blocks[i]:第i块的信息   
void build()  
{  
    int n=100000;  
    block=sqrt(n);  
    num=n/block;  
    if (n%block) num++;  
    for (int i=1;i<=num;i++)  
        l[i]=(i-1)*block+1,r[i]=i*block;  
    r[num]=n;  
    for (int i=1;i<=n;i++)  
        belong[i]=(i-1)/block+1;  
}  
ll ask(int x,int y)  
{  
    ll ret=0;  
    if (belong[x]==belong[y])  
    {  
        for (int i=x;i<=y;i++)  
            ret=(ret+a[i])%mod;  
        return ret;  
    }  
    for (int i=x;i<=r[belong[x]];i++)  
        ret=(ret+a[i])%mod;  
    for (int i=belong[x]+1;i<belong[y];i++)  
        ret=(ret+blocks[i])%mod;  
    for (int i=l[belong[y]];i<=y;i++)  
        ret=(ret+a[i])%mod;  
    return ret;  
}   
```
##### （1）动态区间不同值
```C++
#include <bits/stdc++.h>
#define maxn 500000+5  
#define maxm 1000000+5  
using namespace std;  
int num,belong[maxn],block,l[maxn],r[maxn],last[maxm];  
vector <int> pre[3000];  
int b[maxn],a[maxn];  
//num:分块的个数  
//belong[i]:i所属的块  
//block:块的大小  
//l[i]:i这块的左端点位置  
//r[i]:i这块的右端点位置   
int n,m;   
void build(int n)  
{  
    block=sqrt(n);  
    num=n/block;  
    if (n%block) num++;  
    for (int i=1;i<=num;i++)  
        l[i]=(i-1)*block+1,r[i]=i*block;  
    r[num]=n;  
    for (int i=1;i<=n;i++)  
        belong[i]=(i-1)/block+1;  
}  
void init()  
{  
    scanf("%d%d",&n,&m);  
    build(n);  
    memset(last,0,sizeof(last));  
    for (int i=1;i<=n;i++)  
    {  
        scanf("%d",&a[i]);  
        b[i]=last[a[i]];  
        pre[belong[i]].push_back(b[i]);  
        last[a[i]]=i;  
    }  
    for (int i=1;i<=num;i++)  
        sort(pre[i].begin(),pre[i].end());  
}  
int ask(int x,int y)  
{  
    int ret=0;  
    if (belong[x]==belong[y])  
    {  
        for (int i=x;i<=y;i++)  
            if (b[i]<x) ret++;  
        return ret;  
    }  
    for (int i=x;i<=r[belong[x]];i++)  
        if (b[i]<x) ret++;  
    for (int i=belong[x]+1;i<belong[y];i++)  
        ret+=lower_bound(pre[i].begin(),pre[i].end(),x)-pre[i].begin();  
    for (int i=l[belong[y]];i<=y;i++)  
        if (b[i]<x) ret++;  
    return ret;  
}  
void update(int p,int old,int cur)  
{  
    if (old==cur) return;  
    int pos=0,curnum=belong[p];  
    while (pre[curnum][pos]<old) pos++;  
    pre[curnum][pos]=cur;  
    int size=pre[curnum].size();  
    if (cur>old)  
    {  
        while (pos<size-1 && pre[curnum][pos]>pre[curnum][pos+1])
            swap(pre[curnum][pos],pre[curnum][pos+1]),pos++;  
    }  
    else  
    {  
        while (pos>0 && pre[curnum][pos]<pre[curnum][pos-1])  
            swap(pre[curnum][pos],pre[curnum][pos-1]),pos--;  
    }  
}  
void change(int p,int x)  
{  
    if (a[p]==x) return;  
    int v=a[p],prei;  
    int idx=-1,idv=-1;  
    for (int i=p+1;i<=n;i++)  
    {  
        if (a[i]==x && idx==-1) idx=i;  
        if (a[i]==v && idv==-1) idv=i;  
        if (idx!=-1 && idv!=-1) break;  
    }  
    if (idv!=-1) prei=b[idv],b[idv]=b[p],update(idv,prei,b[idv]);
    if (idx!=-1) prei=b[idx],b[idx]=p,update(idx,prei,b[idx]);  
    idx=-1;  
    for (int i=p-1;i>=1;i--)  
    {  
        if (a[i]==x) idx=i;  
        if (idx!=-1) break;   
    }  
    if (idx!=-1) prei=b[p],b[p]=idx,update(p,prei,b[p]);  
    else prei=b[p],b[p]=0,update(p,prei,b[p]);  
    a[p]=x;  
}  
int main()  
{  
    char s[5];  
    int x,y;  
    init();  
    for (int i=0;i<m;i++)  
    {  
        scanf("%s%d%d",s,&x,&y);  
        x++;  
        if (s[0]=='Q')  
            printf("%d\n",ask(x,y));  
        else change(x,y);  
    }  
    return 0;  
}  
```
#### **6、莫队算法**
```C++
#include <bits/stdc++.h>
using namespace std;
const int N=2e5+5;
int unit,cnt[N],arr[N],res[N],ans=0,post=0;
struct node
{
    int l,r,id;
}q[N];
bool cmp(node a,node b)
{
    return a.l/unit!=b.l/unit?a.l/unit<b.l/unit:a.r<b.r;
}
void add(int pos)
{
    cnt[arr[pos]]++;
    if (arr[pos]>ans && post-cnt[ans]+1>=ans+1)
	{
        post=post-cnt[ans]+1;
        ans++;
    }
    else if (arr[pos]>=ans)
    	post++;
}

void remove(int pos)
{
    cnt[arr[pos]]--;
    if (post==ans && arr[pos]>=ans)
	{
		post=post-1+cnt[ans-1];
        ans--;
    }
    else if (arr[pos]>=ans)
    	post--;
}
void solve(int n,int m)
{
	unit=sqrt(n);
	for (int i=0;i<=n;i++)
		cnt[i]=0;
	ans=post=0;
    for(int i=1;i<=n;i++)
	{	
        scanf("%d",&arr[i]);
    }    
    for(int i=1;i<=m;i++)
	{
        scanf("%d%d",&q[i].l,&q[i].r);
        q[i].id=i;
    }
    sort(q+1,q+m+1,cmp);
    int L=q[1].l,R=L-1;
    for(int i=1;i<=m;i++)
	{
        while(L>q[i].l)
            add(--L);
        while(L<q[i].l)
            remove(L++);
        while(R>q[i].r)
            remove(R--);
        while(R<q[i].r)
            add(++R);
        res[q[i].id]=ans;
    }
    for(int i=1;i<=m;i++)
	{
        printf("%d\n",res[i]);
    }
}
int main()
{
    int n,m;
    while (scanf("%d%d",&n,&m)==2)
    	solve(n,m);
    return 0;
}
```
#### **7、基于哈希的LCP**
```C++
#include <string>  
#include <cstring>  
#include <iostream>  
#include <algorithm>  
#define hb 63  
#define maxl 40000+5  
using namespace std;  
struct Hashnode  
{  
    unsigned long long Hash;  
    int pos;  
}Hash[maxl];  
unsigned long long xp[maxl],h[maxl];  
string s;  
void get_xp()  
{  
    xp[0]=1;  
    for (int i=1;i<maxl;i++)  
        xp[i]=xp[i-1]*hb;  
}  
void get_h()  
{  
    int n=s.length();  
    h[n]=0;  
    h[n-1]=s[n-1];  
    for (int i=n-2;i>=0;i--)  
        h[i]=h[i+1]*hb+s[i];   
}  
int ans,pos,m;  
bool cmp(Hashnode a,Hashnode b)  
{  
    if (a.Hash!=b.Hash) return a.Hash<b.Hash;  
    return a.pos<b.pos;  
}  
bool ok(int L)  
{  
    pos=-1;  
    int n=s.length();  
    for (int i=0;i<n-L+1;i++)  
    {  
        Hash[i].Hash=h[i]-h[i+L]*xp[L];  
        Hash[i].pos=i;  
    }  
    sort(Hash,Hash+n-L+1,cmp);  
    int c;  
    for (int i=0;i<n-L+1;i++)  
    {  
        if ((i==0) || Hash[i].Hash!=Hash[i-1].Hash) c=0;  
        if (++c>=m) pos=max(pos,Hash[i].pos);  
    }  
    return (pos>=0);  
}  
void solve()  
{  
    ans=-1;  
    int l=1,r=s.length();  
    while (l<=r)  
    {  
        int mid=(l+r)/2;  
        if (ok(mid))    l=mid+1;  
        else r=mid-1;  
    }  
    if (r!=0) ans=r;  
    ok(r);   //计算pos位置   
}  
int main()  
{  
    while (cin>>m && m)  
    {  
        cin>>s;  
        get_xp();  
        get_h();  
        solve();  
        if (ans==-1) cout<<"none"<<endl;  
        else cout<<ans<<' '<<pos<<endl;  
    }  
    return 0;  
}  
```
#### **8、二维树状数组**
```C++
//note: 树状数组0是虚拟的节点   
int c[maxn][maxn];  
//note: c[x][y]从[0..x-1][0..y-1]范围内标记和   
int lowbit(int x)  
{  
    return x&(-x);  
}  
void update(int x,int y,int v)  
{  
    for (int i=x;i<maxn;i+=lowbit(i))  
        for (int j=y;j<maxn;j+=lowbit(j))  
            c[i][j]+=v;  
}  
int query(int x,int y)  
{  
    int ret=0;  
    for (int i=x;i>0;i-=lowbit(i))  
        for (int j=y;j>0;j-=lowbit(j))  
            ret+=c[i][j];  
    return ret;  
}  
```
#### **9、最远曼哈顿距离**
最远曼哈顿距离：给点n个mwei坐标 ，求两个坐标$ (x_1,x_2,...,x_m) $ 和$ (y_1,y_2,...,y_m) $
 ，使得这两个点的曼哈顿距离$ |x_1-y_1| + |x_2-y_2| + ... + |x_m-y_m| $最大。
```C++
#include <iostream>  
#include <cstring>  
#include <cstdio>  
#define dim 7    
#define maxn 200000+5  
using namespace std;  
long long p[maxn][dim];  
long long solve(int n)  
{  
    long long minx[1<<dim],maxx[1<<dim];  
    for (int i=0;i<(1<<dim);i++)  
    {  
        minx[i]=(long long)1e15;  
        maxx[i]=(long long)-1e15;  
    }  
    for (int i=0;i<n;i++)  
    {  
        for (int j=0;j<(1<<dim);j++)  
        {  
            long long ret=0;  
            for (int k=0;k<dim;k++)  
            {  
                if (j&(1<<k)) ret+=p[i][k];  
                else ret-=p[i][k];  
            }  
            minx[j]=min(minx[j],ret);  
            maxx[j]=max(maxx[j],ret);  
        }  
    }  
    long long ans=(long long)-1e15;  
    for (int i=0;i<(1<<dim);i++)  
        ans=max(ans,maxx[i]-minx[i]);  
    return ans;  
}  
```
#### **10、带撤销并查集**
```C++
struct DisjointSetUnion   
{
    int fa[maxn*2],rank[maxn*2];
    stack < pair<int*,int> > stk;
    vector <int> history;
    void init()
	{
		for (int i=1;i<=maxn;i++)
			fa[i]=i,rank[i]=0;
	}
    int find(int x)
	{
		return (x^fa[x])?find(fa[x]):x;
	}
    void join(int x,int y)
    {
    		history.push_back(stk.size());
        x=find(x),y=find(y);
        if (x==y) return;
        if (rank[x]<=rank[y])
        {
            stk.push({fa+x,fa[x]}),fa[x]=y;
            if (rank[x]==rank[y]) 
				stk.push({rank+y,rank[y]}),rank[y]++;
        }
        else stk.push({fa+y,fa[y]}),fa[y]=x;
    }
    void undo()
	{
		while (stk.size()!=history.back())
			*stk.top().first=stk.top().second,stk.pop();
		history.pop_back();
	}
    bool check(int x,int y)
	{
		return find(x)==find(y);
	}
}ufs;
```
