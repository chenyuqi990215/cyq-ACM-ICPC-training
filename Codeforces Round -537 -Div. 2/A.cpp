#include <bits/stdc++.h>
using namespace std;
string vowel="aeiou";
string s,t;
int main()
{
	cin>>s>>t;
	if (s.length()!=t.length())
	{
		cout<<"No"<<endl;
		return 0;
	}
	for (int i=0;i<s.length();i++)
	{
		if ((vowel.find(s[i])==string::npos)!=(vowel.find(t[i])==string::npos))
		{
			cout<<"No"<<endl;
			return 0;
		}
	}
	cout<<"Yes"<<endl;
	return 0;
}
