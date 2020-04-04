#include <iostream>
#include <chrono>
#include "polynomial.h"
#include "matrix.h"

using namespace std;
using namespace chrono;
using namespace song;

static const bool close_stream_sync=[]{ios::sync_with_stdio(false);return false;}();
//template class polynomial<complex<double>>;
//template polynomial<complex<double>> from_roots(const vector<complex<double>> &);

int main()
{
    cpolynomial<double> ca={0.0001,0,0,0,0,{-20,-20},{108,12},{-151.5,126.5},
        {-23.75,-134.75},{142.5,-58.75},{-78.75,114.25},{28,-47.75},{-5.5,8.5},1};
    ca=from_roots<complex<double>>({1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2});
//    ca={3,3,2,3,4,1};
    auto t1=high_resolution_clock::now();
    auto x=ca.roots();
    auto t2=high_resolution_clock::now();
//    complex<double> x0=1.0+1e-3;
    complex<double> x0=x[0];
    auto ca1=ca.translation(x0);
    cout<<x0<<"\n";
    for(auto xi:x)
        cout<<xi<<"\n";
//    for(int i=,guard=ca1.size();i<guard;++i)
//        cout<<pow(abs(ca1[0])/abs(ca1[i]),1.0/i)<<"\n";
//    for(int i=0,guard=ca1.size();i<guard;++i)
//        cout<<ca1[i]<<"\n";
    cout<<"calculating time: "<<duration_cast<duration<double,milli>>(t2-t1).count()<<"ms"<<endl;
    return 0;
}
