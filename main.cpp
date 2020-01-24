#include <iostream>
#include <chrono>
#include "polynomial.h"

using namespace std;
using namespace chrono;
using namespace song;

static const bool close_stream_sync=[]{ios::sync_with_stdio(false);return false;}();

int main()
{
    cpolynomial<double> ca={0.0001,0,0,0,0,{-20,-20},{108,12},{-151.5,126.5},
        {-23.75,-134.75},{142.5,-58.75},{-78.75,114.25},{28,-47.75},{-5.5,8.5},1};
//    ca=from_roots<complex<double>>({1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2});
    auto t1=high_resolution_clock::now();
    auto x=ca.roots();
    auto t2=high_resolution_clock::now();
    for(const auto &c:x)
        cout<<c<<" "<<std::abs(ca.offset(c))<<" "<<std::abs(ca(c))<<"\n";
    ca=from_roots(x);
    for(auto &c:ca)
        cout<<c<<" ";
//    cout<<x<<" "<<ca(x)<<endl;
    cout<<"calculating time: "<<duration_cast<nanoseconds>(t2-t1).count()/1.e6<<"ms"<<endl;
    return 0;
}
