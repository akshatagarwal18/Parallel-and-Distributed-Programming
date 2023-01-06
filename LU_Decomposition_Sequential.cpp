#include<bits/stdc++.h>
using namespace std;
using ll = long long ;
using ld = long double;
int main(){
ll n; cin>>n;
ld a[n][n];
ld l[n][n];
ld u[n][n];
ld res[n][n];
for(int i = 0 ; i < n ; i++){
    for(int j = 0 ; j < n ; j++){
        cin>>a[i][j];
        if(i==j){l[i][j]=1.0;}
        if(i>j){u[i][j]=0.0;}
        if(i<j){l[i][j]=0.0;}
        res[i][j]=0.0;
    }
}

ll p[n]; for(int i = 0 ; i < n ; i++){p[i]=i;}

for(int k = 0 ; k < n ; k++){
    ld mx = 0.0;
    ll k1 = -1;
    for(int i = k ; i < n ; i++){
        ld temp = abs(a[i][k]);
        if(temp > mx){
            mx=temp;
            k1=i;
        }
    }

    swap(p[k],p[k1]);
    for(int i = 0 ; i < n ; i++){
        swap(a[k][i],a[k1][i]);
    }
    for(int i = 0 ; i < k ; i++){
        swap(l[k][i],l[k1][i]);
    }
    u[k][k]=a[k][k];
    for(int i = k+1; i < n ; i++){
        l[i][k] = a[i][k]/u[k][k];
        u[k][i] = a[k][i];
    }

    for(int i = k+1; i < n ; i++){
        for(int j = k+1; j < n ; j++){
            a[i][j] = a[i][j] - l[i][k]*u[k][j]; 
        }
    }

}

for(int i = 0 ; i < n ; i++){cout<<p[i]<<" ";}
cout<<endl<<endl;
for(int i = 0 ; i < n ; i++){
    for(int j = 0 ; j < n ; j++){
        cout<<l[i][j]<<" ";
    }cout<<endl;
}
cout<<endl<<endl;
for(int i = 0 ; i < n ; i++){
    for(int j = 0 ; j < n ; j++){
        cout<<u[i][j]<<" ";
    }cout<<endl;
}
cout<<endl<<endl;
for(int i = 0 ; i < n ; i++){
    for(int j = 0 ; j < n ; j++){
        for(int k = 0 ; k < n ; k++){
            res[i][j]+=(l[i][k]*u[k][j]);
        }
    }
}
for(int i = 0 ; i < n ; i++){
    for(int j = 0 ; j < n ; j++){
        cout<<res[i][j]<<" ";
    }cout<<endl;
}


return 0;    
}
