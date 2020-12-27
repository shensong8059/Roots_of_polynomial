#include <iostream>
#include <chrono>
#include "polynomial.h"
#include "matrix.h"

using namespace std;
using namespace chrono;
using namespace song;

namespace
{
    const bool close_stream_sync=[]{ios::sync_with_stdio(false);return false;}();
}
//template class polynomial<complex<double>>;
//template polynomial<complex<double>> from_roots(const vector<complex<double>> &);

int main()
{
    matrix<float> m={3,3,{5,-3,2,6,-4,4,4,-4,5}};
    m=m.Hessenberg();
    m.HessenQR();
//    cpolynomial<float> ca={0.0001,0,0,0,0,{-20,-20},{108,12},{-151.5,126.5},
//        {-23.75,-134.75},{142.5,-58.75},{-78.75,114.25},{28,-47.75},{-5.5,8.5},1};
    auto ca=from_roots<complex<float>>({0,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2});
//    ca={3,3,2,3,4,1};
    auto t1=high_resolution_clock::now();
//    polynomial<double> cb{1,3,3,1};
//    auto vb=cb.roots_of_degree3();
    auto x=ca.roots();
    auto t2=high_resolution_clock::now();
    for(auto xi:x)
        cout<<xi<<'\t'<<std::abs(ca(xi))<<"\n";
    cout<<"calculating time: "<<duration_cast<duration<double,milli>>(t2-t1).count()<<"ms"<<endl;
    return 0;
}
