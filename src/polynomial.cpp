#include "polynomial.h"

namespace song
{
    using namespace std;
    template<class T>
    vector<typename polynomial<T>::root_type> polynomial<T>::roots()const
    requires polynomial<T>::complex_floating_point
    {
//            static_assert(is_std_complex_v<coefficient_type>,
//                          "Can not get roots of non-complex coefficient polynomial");
        int n=degree();
        if(n<=0)
            throw std::runtime_error("degree is less than zero");
        auto an=this->back();
        if(std::abs(an)<eps)
            throw std::runtime_error("first term is too small");
        constexpr std::vector<coefficient_type> (polynomial::*reg_rt_memfunc[])()const=
        {
            &polynomial::roots_of_degree1,
            &polynomial::roots_of_degree2,
            &polynomial::roots_of_degree3,
            &polynomial::roots_of_degree4
        };
        int rs=std::size(reg_rt_memfunc);
        if(n<=rs)
            return (this->*reg_rt_memfunc[n-1])();
        auto f=*this;
        auto inv_an=base_type(1.0)/this->back();
        for(auto &c:f)
            c*=inv_an;
        std::vector<coefficient_type> ans;
        auto f0=f;
        ans.reserve(n);
        while(f.degree()>rs)
        {
            auto temp=f.root_without_init();
            ans.push_back(temp);
            f/={-temp,1};
        }
        auto last_ans=(f.*reg_rt_memfunc[rs-1])();
        for(auto &lans:last_ans)
            ans.push_back(std::move(lans));
        for(auto &x:ans)
            x=f0.root_with_init(x);
        return ans;
    }

    template vector<polynomial<complex<float>>::root_type> polynomial<complex<float>>::roots()const;
    template vector<polynomial<complex<double>>::root_type> polynomial<complex<double>>::roots()const;
    template vector<polynomial<complex<long double>>::root_type> polynomial<complex<long double>>::roots()const;
//    template vector<polynomial<complex<double>>::root_type> polynomial<complex<double>>::roots()const;
//    template vector<polynomial<complex<double>>::root_type> polynomial<complex<double>>::roots()const;
//    template vector<polynomial<complex<double>>::root_type> polynomial<complex<double>>::roots()const;
}
