#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <iterator>//back_inserter_iterator
#include <type_traits>
#include <concepts>
#include <numbers>

namespace song
{
    template<class T>
    requires std::floating_point<T>||
    (std::same_as<T,std::complex<typename T::value_type>>&&
        std::floating_point<typename T::value_type>)
    class polynomial:public std::vector<T>
    {
        template<class U>
        struct is_std_complex:std::false_type{};
        template<class U>
        struct is_std_complex<std::complex<U>>:std::true_type{};
        template<class U>
        constexpr static bool is_std_complex_v=is_std_complex<U>::value;
    public:
        typedef T coefficient_type;
        using abs_coefficient_type=decltype(std::abs(coefficient_type()));
        static constexpr auto eps=std::numeric_limits<abs_coefficient_type>::epsilon();//epsilon()==2.22045e-16;
        static constexpr auto inf=std::numeric_limits<abs_coefficient_type>::infinity();
        static constexpr auto PI=std::numbers::pi_v<abs_coefficient_type>;
    public:
        coefficient_type root_with_init(const coefficient_type &arg_x)const
        {
            auto x=arg_x;
            auto dx=offset(x);
            auto err=std::abs((*this)(x));
            while(err>eps&&std::abs(dx)>eps)
            {
                auto cur_x=x+dx;
                auto err1=std::abs((*this)(cur_x));
                if(err1>err)
                {
                    break;
                }
                x=cur_x;
                dx=offset(x);
                err=err1;
            }
            return x;
        }
        coefficient_type root_without_init()const
        {
            int deg=degree();
            if(deg<=0)
                throw std::runtime_error("degree is less than zero");
            coefficient_type a=this->front()/this->back();
            abs_coefficient_type h=1.0/deg;
            abs_coefficient_type x_r=std::pow(std::abs(a),h);
            coefficient_type x=0,dx=offset(x);
            abs_coefficient_type err=std::abs((*this)(dx));
            for(int i=0;i<deg;++i)
            {
                auto curx=std::polar(x_r,2*i*(PI*h));
                auto curdx=offset(curx);
                auto curerr=std::abs((*this)(curdx));
                if(curerr<err)
                {
                    x=curx;
                    dx=curdx;
                    err=curerr;
                }
            }
            return root_with_init(x+dx);
        }
//    private:
        coefficient_type operator()(const coefficient_type &x)const
        {
            if(this->empty())
                return coefficient_type(0.0);
            int n=this->degree();
            coefficient_type p=(*this)[n];
            for(int j=n-1;j>=0;--j)
                p=p*x+(*this)[j];
            return p;
        }
        coefficient_type div_on_point(const polynomial &g,const coefficient_type &z)const
        {
            if(g.empty())
                throw std::runtime_error("g is zero in polynomial::div_on_point");
            auto fz=(*this)(z),gz=g(z);
            if(std::abs(gz)>eps)
                return fz/gz;
            if(std::abs(fz)>eps)
                return {inf,inf};
            auto [q,r]=poly_divide(g);
            auto dr=r.derivate(),dg=g.derivate();
            while(true)
            {
                auto dgz=dg(z),drz=dr(z);
                if(std::abs(dgz)>eps)
                {
                    return q(z)+drz/dgz;
                }
                if(std::abs(drz)>eps)
                {
                    return {inf,inf};
                }
                auto t=dg.back();
                dr/=t;dg/=t;
                dr=dr.derivate();
                dg=dg.derivate();
            }
            return {};
        }
        polynomial translation(const coefficient_type &x)const
        {
            auto dfi_div_ifact=*this;
            polynomial ret(this->size());
            ret[0]=dfi_div_ifact(x);
            for(int i=1,guard=ret.size();i<guard;++i)
            {
                dfi_div_ifact=dfi_div_ifact.derivate();
                dfi_div_ifact/=coefficient_type(i);
                ret[i]=dfi_div_ifact(x);
            }
            return ret;
        }
        std::vector<coefficient_type> roots_of_degree1()const
        {
            return {-(*this)[0]/(*this)[1]};
        }
        std::vector<coefficient_type> roots_of_degree2()const
        {
            auto a=(*this)[2],b=(*this)[1],c=(*this)[0];
            auto p=b/a,q=c/a;
            auto sqrt_delta=std::sqrt(p*p-4.0*q);
            return {(-p+sqrt_delta)/2.0,(-p-sqrt_delta)/2.0};
        }
        //from https://zhuanlan.zhihu.com/p/40349993
        std::vector<coefficient_type> roots_of_degree3()const
        {
            auto a=(*this)[3];
            auto one_div_a=1.0/a;
            auto b=(*this)[2]*one_div_a,c=(*this)[1]*one_div_a,d=(*this)[0]*one_div_a;
            auto b_2=b*b;
            auto p=c-1.0/3*b_2;
            auto q=d-1.0/3*b*c+2.0/27.0*b_2*b;
            auto w=std::polar(1.,2/3.0*PI);
            auto w2=std::polar(1.,-2/3.0*PI);
            auto q1=q/2.0,p1=1.0/3*p,b1=1.0/3*b;
            auto sq_delt=std::sqrt(q1*q1+p1*p1*p1);
            auto t1=std::pow(-q1+sq_delt,1.0/3),t2=std::pow(-q1-sq_delt,1.0/3);
            return {t1+t2-b1,w*t1+w2*t2-b1,w2*t1+w*t2-b1};
        }
        //from https://www.cnblogs.com/larissa-0464/p/11706131.html
        std::vector<coefficient_type> roots_of_degree4()const
        {
            auto a=(*this)[4];
            auto one_div_a=1.0/a;
            auto b=(*this)[3]*one_div_a;
            auto c=(*this)[2]*one_div_a;
            auto d=(*this)[1]*one_div_a;
            auto e=(*this)[0]*one_div_a;
            auto P=(c*c+12.0*e-3.0*b*d)/9.0;
            auto Q = (27.0*d*d+2.0*c*c*c+27.*b*b*e-72.*c*e-9.*b*c*d)/54.;
            auto D = std::sqrt(Q*Q-P*P*P);
            auto t1=Q+D,t2=Q-D;
            auto u=std::abs(t1)>std::abs(t2)?std::pow(t1,1.0/3):std::pow(t2,1.0/3);
            auto v=std::sqrt(std::abs(u))<eps?coefficient_type(0.):P/u;
            coefficient_type w[]={std::polar(1.,2./3*PI),std::polar(1.,-2./3*PI)};
            auto temp0=b*b-8./3*c;
            auto temp1=4.*(u+v);
            coefficient_type m2=temp0+temp1;
            for(int i=0;i<2;++i)
            {
                auto tmp1=4.*(w[i]*u+w[1-i]*v);
                auto temp=temp0+tmp1;
                if(std::abs(temp)<std::abs(temp))
                {
                    m2=temp;
                    temp1=tmp1;
                }
            }
            auto m=std::sqrt(m2);
            coefficient_type S,TT;
            if(std::abs(m)<eps)
            {
                S=temp0;
                TT=0.;
            }
            else
            {
                S=2.*temp0-temp1;
                TT=(8.*b*c-16.*d-2.*b*b*b)/m;
            }
            auto bm0=-b-m,bm1=-b+m;
            auto sq_ST0=std::sqrt(S-TT),sq_ST1=std::sqrt(S+TT);
            return {(bm0+sq_ST0)/4.,(bm0-sq_ST0)/4.,(bm1+sq_ST1)/4.,(bm1-sq_ST1)/4.};
        }
        coefficient_type offset(const coefficient_type &x)const
        {
            auto dpn=derivate(),d2pn=dpn.derivate();
            auto pnx=(*this)(x),dpnx=dpn(x),d2pnx=d2pn(x);
            int n=degree();
            auto flag1=std::abs(dpnx)<eps,flag0=std::abs(pnx)<eps;
            if(!flag1)
            {
                auto t=1.0/dpnx,temp1=pnx*t,temp2=temp1*d2pnx*t;
                auto delta=std::sqrt(abs_coefficient_type(n-1)*(abs_coefficient_type(n-1)-abs_coefficient_type(n)*temp2));
                auto d1=1.0+delta,d2=1.0-delta;
                return -temp1*(abs_coefficient_type(n)/(std::abs(d1)>std::abs(d2)?d1:d2));
            }
            // flag1
            if(!flag0)
            {
                //dpnx==0,pnx!=0,
                if(std::abs(d2pnx)<eps)
                    return {inf,inf};
                auto delta=std::sqrt(abs_coefficient_type(n-1)*(abs_coefficient_type(n-1)*(dpnx*dpnx)-abs_coefficient_type(n)*(pnx*d2pnx)));
                auto d1=dpnx+delta,d2=dpnx-delta;
                return -abs_coefficient_type(n)*(pnx/(std::abs(d1)>std::abs(d2)?d1:d2));
            }
            // flag1&&flag2
            auto temp1=div_on_point(dpn,x);
            auto temp2=poly_multiply(d2pn).div_on_point(dpn.poly_multiply(dpn),x);
            return -temp1/(1.0-temp2);
        }

    public:
        using std::vector<T>::vector;//继承vector构造函数，不包括从vector到poly的转换
        polynomial(const std::vector<coefficient_type> &p):std::vector<coefficient_type>(p){}//实例化要补全类型参数
        polynomial(std::vector<coefficient_type> &&p):std::vector<coefficient_type>(std::move(p)){}

        int degree()const
        //最高次数，值为-1表示零多项式
        {
            return int(this->size())-1;
        }
        std::vector<coefficient_type> roots()const
        requires is_std_complex_v<coefficient_type>
        {
//            static_assert(is_std_complex_v<coefficient_type>,
//                          "Can not get roots of non-complex coefficient polynomial");
            int n=degree();
            if(n<=0)
                throw std::runtime_error("degree is less than zero");
            auto a0=this->back();
            if(std::pow(std::abs(a0),1./n)<eps)
                throw std::runtime_error("first term is too small");
            constexpr std::vector<coefficient_type> (polynomial::*reg_rt_memfunc[])()const=
            {
                &polynomial::roots_of_degree1,
                &polynomial::roots_of_degree2,
                &polynomial::roots_of_degree3,
                &polynomial::roots_of_degree4
            };
            if(n<=4)
                return (this->*reg_rt_memfunc[n-1])();
            auto f=*this;
            std::vector<coefficient_type> ans;
            ans.reserve(n);
            while(f.degree()>4)
            {
                auto temp=f.root_without_init();
                ans.push_back(temp);
                f=f.div_monomial_factor(temp).first;
            }
            auto last_ans=f.roots_of_degree4();
            move(last_ans.cbegin(),last_ans.cend(),std::back_insert_iterator(ans));
            for(auto &x:ans)
                x=root_with_init(x);
            return ans;
        }
        polynomial mul_monomial_factor(const coefficient_type &a)const
        //f*(x-a)
        {
            int n=this->size();
            polynomial c(n+1);
            c[n]=(*this)[n-1];
            for(int j=n-1;j>=1;--j)
            {
                c[j]=(*this)[j-1]-(*this)[j]*a;
            }
            c[0]=(*this)[0]*(-a);
            return c;
        }
        polynomial &operator%=(const polynomial &v)
        {
            this->swap(poly_divide(v).second);
            return *this;
        }
        polynomial &operator*=(const polynomial &rhs)
        {
            this->swap(poly_multiply(rhs));
            return *this;
        }
        polynomial &operator+=(const polynomial &rhs)
        {
            const int m=this->size(),n=rhs.size();
            if(m<n)
            {
                for(int i=0;i<m;++i)
                {
                    (*this)[i]+=rhs[i];
                }
                this->reserve(n);
                for(int i=m;i<n;++i)
                {
                    this->push_back(rhs[i]);
                }
            }
            else
            {
                for(int i=0;i<n;++i)
                {
                    (*this)[i]+=rhs[i];
                }
            }
            return normalize();
        }
        polynomial &operator-=(const polynomial &rhs)
        {
            const int m=this->size(),n=rhs.size();
            if(m<n)
            {
                for(int i=0;i<m;++i)
                {
                    (*this)[i]-=rhs[i];
                }
                this->reserve(n);
                for(int i=m;i<n;++i)
                {
                    this->push_back(rhs[i]);
                }
            }
            else
            {
                for(int i=0;i<n;++i)
                {
                    (*this)[i]-=rhs[i];
                }
            }
            return normalize();
        }
        polynomial &operator/=(const polynomial &v)
        {
            this->swap(poly_divide(v).first);
            return *this;
        }
        polynomial &operator*=(const coefficient_type &c)
        {
            for(auto &pi:*this)
                pi*=c;
            monic();
            return *this;
        }
        polynomial &operator/=(const coefficient_type &c)
        {
            return (*this)*=1.0/c;
        }

        polynomial operator-()&&
        {
            for(auto &c:*this)
            {
                c=-c;
            }
            return std::move(*this);
        }
        polynomial operator-()const&
        {
            return -(polynomial(*this));
        }
        //规范化，不允许0系数占据最高次
        polynomial &normalize()
        {
            this->erase(std::find_if(this->rbegin(),this->rend(),[](const coefficient_type &x){return std::abs(x)>eps;}).base(),this->end());
            return *this;
        }
        polynomial monic()const
        {
            polynomial ret(this->size());
            auto one_d_an=coefficient_type(1.0)/this->back();
            for(int i=0,guard=this->size();i<guard;++i)
            {
                ret[i]=(*this)[i]*one_d_an;
            }
            return ret;
        }
        std::pair<polynomial,polynomial> poly_divide(const polynomial &v)const
        {
            if(v.empty())
                throw std::runtime_error("Divisor is close to zero");
            int n=this->size()-1,nv=v.size()-1;
            polynomial r=*this,q(n-nv+1);
            for(int k=n-nv;k>=0;--k)
            {
                q[k]=r[nv+k]/v[nv];
                for(int j=nv+k-1;j>=k;--j)
                    r[j]-=q[k]*v[j-k];
            }
            r.erase(r.cbegin()+nv,r.cend());
            r.normalize();
            return std::pair{std::move(q),std::move(r)};
        }
        polynomial poly_multiply(const polynomial &rhs)const
        {
            int m=this->size(),n=rhs.size();
            if(m<n)
                return rhs.poly_multiply(*this);
            int p=m+n-1;
            polynomial c(p);
            for(int i=0;i<n;++i)
            {
                coefficient_type ci(0.0);
                for(int j=0;j<=i;++j)
                    ci+=rhs[j]*(*this)[i-j];
                c[i]=ci;
            }
            for(int i=n;i<m;++i)
            {
                coefficient_type ci(0.0);
                for(int j=0;j<n;++j)
                    ci+=rhs[j]*(*this)[i-j];
                c[i]=ci;
            }
            for(int i=m;i<p;++i)
            {
                coefficient_type ci(0.0);
                for(int j=i-m;j<n;++j)
                    ci+=rhs[j]*(*this)[i-j];
                c[i]=ci;
            }
            return c;
        }
        polynomial poly_add(const polynomial &rhs)const
        {
            auto p=*this;
            p+=rhs;
            return p;
        }
        polynomial poly_substract(const polynomial &rhs)const
        {
            auto p=*this;
            p-=rhs;
            return p;
        }
        polynomial derivate()const
        {
            int n=this->size();
            if(n<=1)
                return {};
            polynomial p(n-1);
            for(int i=1;i<n;++i)
            {
                p[i-1]=abs_coefficient_type(i)*(*this)[i];
            }
            return p;
        }

        std::pair<polynomial,coefficient_type> div_monomial_factor(const coefficient_type &a)const
        {
            int n=degree();
            coefficient_type rem=(*this)[n];
            polynomial quotient(n);
            for(int i=n-1;i>=0;--i)
            {
                quotient[i]=rem;
                rem=(*this)[i]+rem*a;
            }
            return std::pair{std::move(quotient),rem};
        }
    };

    template<class T>
    inline polynomial<T> operator+(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.polyadd(rhs);
    }

    template<class T>
    inline polynomial<T> operator-(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.polysub(rhs);
    }

    template<class T>
    inline polynomial<T> operator/(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.poly_divide(rhs).first;
    }

    template<class T>
    inline polynomial<T> operator%(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.poly_divide(rhs).second;
    }

    template<class T>
    inline polynomial<T> operator*(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.poly_multiply(rhs);
    }

    template<class T>
    inline polynomial<T> from_roots(const std::vector<T> &roots)
    {
        polynomial<T> ret={1.0};
        for(auto &c:roots)
        {
            ret=ret.mul_monomial_factor(c);
        }
        return ret;
    }

    template<class T>
    using cpolynomial=polynomial<std::complex<T>>;
}

#endif // POLYNOMIAL_H
