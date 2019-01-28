#include <bits/stdc++.h>
using namespace std;
int main()
{
	string s;
	int k;
	cin>>s>>k;
	int lmin=0,lmax=0,p=0;
	for (int i=0;i<s.length();i++)
	{
		if (isalpha(s[i])) lmin++,lmax++,p++;
		if (s[i]=='?') lmin--;
		if (s[i]=='*') lmin--,lmax=200;
	}
	if (k<lmin || k>lmax) printf("Impossible\n");
	else
	{
		if (k<=p)
		{
			string ans="";
			int t=p-k;
			for (int i=0;i<s.length();i++)
			{
				if (s[i]=='?' || s[i]=='*')
				{
					if (t>0)
					{
						t--;
						ans=ans.substr(0,ans.length()-1);
					}
				}
				else ans+=s[i];
			}
			cout<<ans<<endl;
		}
		else
		{
			string ans="";
			int t=k-p;
			bool flag=true;
			for (int i=0;i<s.length();i++)
			{
				if (s[i]=='*')
				{
					if (flag)
					{
						for (int j=0;j<t;j++)
							ans+=s[i-1];
						flag=false;	
					}
				}
				else if (s[i]!='?')
					ans+=s[i];
			}
			cout<<ans<<endl;
		}
	}
	return 0;
}
