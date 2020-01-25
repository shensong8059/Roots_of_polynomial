#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cmath>
#include <complex>

namespace song
{
    template<class T>
    class basic_polynomial:public std::vector<T>
    {
        template<class U>
        struct base
        {
            typedef U type;
        };
        template<class U>
        struct base<std::complex<U>>
        {
            typedef U type;
        };
    public:
        typedef T coef_type;
        typedef typename base<coef_type>::type base_type;
        static constexpr base_type eps=1000*std::numeric_limits<base_type>::epsilon();//2.22045e-13;
        static constexpr auto inf=std::numeric_limits<base_type>::infinity();
    protected:
        using std::vector<coef_type>::vector;//继承vector构造函数，不包括从vector到poly的转换
        template<template<class> class Container>
        static auto fix_iterator(const Container<coef_type> &v)
        {
            auto it_b=v.begin(),it_e=std::prev(v.end());
            auto i=std::distance(it_b,it_e);
            for(;it_e!=it_b;--it_e)
            {
                auto x=std::pow(std::abs(*it_e),1.0/i);
                if(x>eps)
                    return std::pair{it_b,std::next(it_e)};
                --i;
            }

            if(std::abs(*it_b)>eps)
                ++it_e;
            return std::pair{it_b,it_e};
        }
        template<template<class> class Container>
        static auto fix_container(const Container<coef_type> &v)
        {
            auto temp=fix_iterator(v);
            return std::vector<coef_type>{temp.first,temp.second};
        }
        static auto fix_container(std::vector<coef_type> &&v)
        {
            auto temp=fix_iterator(v);
            v.erase(temp.second,v.end());
            return v;
        }
        basic_polynomial &correct()
        {
            auto temp=fix_iterator(*this);
            this->erase(temp.second,this->end());
            return *this;
        }
        basic_polynomial translation(const coef_type &x)const
        {
            auto dfi_div_ifact=*this;
            std::vector<coef_type> ret(this->size());
            ret[0]=dfi_div_ifact(x);
            for(int i=1,guard=ret.size();i<guard;++i)
            {
                dfi_div_ifact=dfi_div_ifact.derivate();
                dfi_div_ifact/=coef_type(i);
                ret[i]=dfi_div_ifact(x);
            }
            return ret;
        }
    public:
        basic_polynomial(const std::vector<coef_type> &v):std::vector<coef_type>(fix_container(v)){}
        basic_polynomial(std::vector<coef_type> &&v):std::vector<coef_type>(fix_container(std::move(v))){}
        basic_polynomial(std::initializer_list<coef_type> il):std::vector<coef_type>(fix_container(il)){}
        coef_type operator()(const coef_type &x)const
        {
            if(this->empty())
                return base_type(0.0);
            const int n=this->size()-1;
            coef_type p=(*this)[n];
            for(int j=n-1;j>=0;--j)
                p=p*x+(*this)[j];
            return p;
        }
        basic_polynomial &operator+=(const basic_polynomial &rhs)
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
            return this->correct();
        }
        basic_polynomial &operator-=(const basic_polynomial &rhs)
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
            return this->correct();
        }
        basic_polynomial mul_monomial_factor(const coef_type &a)const
        //f*(x-a)
        {
            int n=this->size();
            basic_polynomial c(n+1);
            c[n]=(*this)[n-1];
            for(int j=n-1;j>=1;--j)
            {
                c[j]=(*this)[j-1]-(*this)[j]*a;
            }
            c[0]=(*this)[0]*(-a);
            return c;
        }
        auto div_monomial_factor(const coef_type &a)const
        //f/(x-a)
        {
            int n=this->size()-1;
            coef_type rem=(*this)[n];
            basic_polynomial quotient(n);
            for(int i=n-1;i>=0;--i)
            {
                quotient[i]=rem;
                rem=(*this)[i]+rem*a;
            }
            return std::pair{quotient,rem};
        }
        auto polydiv(const basic_polynomial &v)const
        //f/g
        {
            int n=this->size()-1,nv=v.size()-1;
            basic_polynomial r=*this,q(n);
            for(int k=n-nv;k>=0;--k)
            {
                q[k]=r[nv+k]/v[nv];
                for(int j=nv+k-1;j>=k;--j)
                    r[j]-=q[k]*v[j-k];
            }
            r.erase(r.cbegin()+nv,r.cend());
            r.correct();
            return std::pair{q,r};
        }
        basic_polynomial &operator/=(const basic_polynomial &v)
        {
            this->swap(this->polydiv(v).first);
            return *this;
        }
        basic_polynomial &operator%=(const basic_polynomial &v)
        {
            this->swap(this->polydiv(v).second);
            return *this;
        }
        basic_polynomial polymul(const basic_polynomial &rhs)const
        //f*g
        {
            int m=this->size(),n=rhs.size();
            if(m<n)
                return rhs.polymul(*this);
            int p=m+n-1;
            basic_polynomial c(p);
            for(int i=0;i<n;++i)
            {
                coef_type ci(0.0);
                for(int j=0;j<=i;++j)
                    ci+=rhs[j]*(*this)[i-j];
                c[i]=ci;
            }
            for(int i=n;i<m;++i)
            {
                coef_type ci(0.0);
                for(int j=0;j<n;++j)
                    ci+=rhs[j]*(*this)[i-j];
                c[i]=ci;
            }
            for(int i=m;i<p;++i)
            {
                coef_type ci(0.0);
                for(int j=i-m;j<n;++j)
                    ci+=rhs[j]*(*this)[i-j];
                c[i]=ci;
            }
            return c;
        }
        basic_polynomial polyadd(const basic_polynomial &rhs)const
        //f+g
        {
            auto p=*this;
            return p+=rhs;
        }
        basic_polynomial polysub(const basic_polynomial &rhs)const
        //f-g
        {
            auto p=*this;
            return p-=rhs;
        }
        basic_polynomial derivate()const
        //df/dx
        {
            int n=this->size();
            if(n<=1)
                return {};
            basic_polynomial p(n-1);
            for(int i=1;i<n;++i)
            {
                p[i-1]=base_type(i)*(*this)[i];
            }
            return p;
        }
        basic_polynomial &operator*=(const basic_polynomial &rhs)
        {
            return *this=this->polymul(rhs);
        }
        basic_polynomial &operator*=(const coef_type &c)
        {
            for(auto &pi:*this)
                pi*=c;
            return *this;
        }
        basic_polynomial &operator/=(const coef_type &c)
        {
            return (*this)*=1.0/c;
        }
        int degree()const
        //最高次数，值为-1表示零多项式
        {
            return int(this->size())-1;
        }
    };
    template<template<class T> class Container,class T>
    inline basic_polynomial<T> from_roots(const Container<T> &roots)
    {
        basic_polynomial<T> ret={1.0};
        for(auto &c:roots)
        {
            ret=ret.mul_monomial_factor(c);
        }
        return ret;
    }
    template<class T>
    inline basic_polynomial<T> from_roots(std::initializer_list<T> roots)
    {
        basic_polynomial<T> ret={1.0};
        for(auto &c:roots)
        {
            ret=ret.mul_monomial_factor(c);
        }
        return ret;
    }
    template<class T>
    class polynomial:public basic_polynomial<T>
    {
    public:
        using basic_polynomial<T>::basic_polynomial;
        polynomial(const basic_polynomial<T> &rhs):basic_polynomial<T>(rhs){}
        polynomial(basic_polynomial<T> &&rhs):basic_polynomial<T>(std::move(rhs)){}
    };

    template<class T>
    class cpolynomial:public basic_polynomial<std::complex<T>>
    {
    public:
        typedef typename basic_polynomial<std::complex<T>>::coef_type coef_type;
        typedef typename basic_polynomial<coef_type>::base_type base_type;
        static constexpr base_type PI=4*std::atan(base_type(1.0));
        static constexpr auto eps=basic_polynomial<std::complex<T>>::eps;
    private:
        coef_type root_on_guess(const coef_type &arg_x)const
        {
            auto x=arg_x;
            auto dx=this->offset(x);
            auto err=std::abs((*this)(x));
            while(std::abs(dx)>this->eps)
            {
                auto cur_x=x-dx;
                auto err1=std::abs((*this)(cur_x));
                if(err1>err)
                {
                    break;
                }
                x=cur_x;
                dx=this->offset(x);
                err=err1;
            }
            dx=this->tiny_offset(x);
            err=std::abs((*this)(x));
            while(std::abs(dx)>this->eps)
            {
                auto cur_x=x+dx;
                auto err1=std::abs((*this)(cur_x));
                if(err1>err)
                {
                    break;
                }
                x=cur_x;
                dx=this->tiny_offset(x);
                err=err1;
            }
            return x;
        }
        coef_type root_without_init()const
        {
            int deg=this->degree();
            if(deg<=0)
                throw std::runtime_error("degree is less than zero");
            coef_type a=this->front()/this->back();
            base_type h=1.0/deg;
            base_type x_r=std::pow(std::abs(a),h);
            coef_type x=0,dx=this->offset(x);
            base_type err=std::abs((*this)(dx));
            for(int i=0;i<deg;++i)
            {
                auto curx=std::polar(x_r,2*i*(PI*h));
                auto curdx=this->offset(curx);
                auto curerr=std::abs((*this)(curdx));
                if(curerr<err)
                {
                    x=curx;
                    dx=curdx;
                    err=curerr;
                }
            }
            return root_on_guess(x-dx);
        }
        coef_type div_on_point(const cpolynomial &g,const coef_type &z)const
        {
            if(g.empty())
                throw std::runtime_error("g is zero in cpolynomial<...>::div_on_point");
            auto fz=(*this)(z),gz=g(z);
            if(std::abs(gz)>eps)
                return fz/gz;
            if(std::abs(fz)>eps)
                return {this->inf,this->inf};
            auto [q,r]=this->polydiv(g);
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
                    return {this->inf,this->inf};
                }
                auto t=dg.back();
                dr/=t;dg/=t;
                dr=dr.derivate();
                dg=dg.derivate();
            }
            return {};
        }
        cpolynomial translation(const coef_type &x)const
        {
            return this->basic_polynomial<coef_type>::translation(x);
        }
    public:
        using basic_polynomial<coef_type>::basic_polynomial;
        cpolynomial(const basic_polynomial<coef_type> &p):basic_polynomial<coef_type>(p){}//实例化要补全类型参数
        cpolynomial(basic_polynomial<coef_type> &&p):basic_polynomial<coef_type>(std::move(p)){}
        coef_type offset(const coef_type &x)const
        {
            auto dpn=this->derivate(),d2pn=dpn.derivate();
            auto pnx=(*this)(x),dpnx=dpn(x),d2pnx=d2pn(x);
            int n=this->degree();
            auto flag1=std::abs(dpnx)<eps,flag0=std::abs(pnx)<eps;
            if(!flag1)
            {
                auto t=1.0/dpnx,temp1=pnx*t,temp2=temp1*d2pnx*t;
                auto delta=std::sqrt(base_type(n-1)*(base_type(n-1)-base_type(n)*temp2));
                auto d1=1.0+delta,d2=1.0-delta;
                return temp1*(base_type(n)/(std::abs(d1)>std::abs(d2)?d1:d2));
            }
            // flag1
            if(!flag0)
            {
                if(std::abs(d2pnx)<eps)
                    return {this->inf,this->inf};
                auto delta=std::sqrt(base_type(n-1)*(base_type(n-1)*(dpnx*dpnx)-base_type(n)*(pnx*d2pnx)));
                auto d1=dpnx+delta,d2=dpnx-delta;
                return base_type(n)*(pnx/(std::abs(d1)>std::abs(d2)?d1:d2));
            }
            // flag1&&flag2
            auto temp1=this->div_on_point(dpn,x);
            auto temp2=this->polymul(d2pn).div_on_point(dpn.polymul(dpn),x);
            return temp1/(1.0-temp2);
        }
        coef_type tiny_offset(const coef_type &x)const
        {
            auto f=translation(x);
//            std::cout<<x<<std::endl;
//            for(auto &c:f)
//                std::cout<<std::abs(c)<<" "<<std::flush;
//            std::cout<<std::endl;

            auto f0=f.front();
            auto af0=std::abs(f0);
            auto xx=this->inf;
            int k;
            for(int i=1,guard=f.size();i<guard;++i)
            {
                auto afi=std::abs(f[i]);
                if(std::pow(afi,1.0/i)>this->eps)
                {
                    auto xi=std::pow(af0/afi,1.0/i);
                    if(xi<xx)
                    {
                        xx=xi;
                        k=i;
                    }
                }
            }
//            std::cout<<k<<std::endl;

            auto cxk=std::pow(-f0/f[k],1.0/k);
            auto d=std::abs(cxk),alpha=std::arg(cxk);
            auto err=this->inf;
            coef_type ret;
            for(int i=0;i<k;++i)
            {
                auto c=std::polar(d,alpha+2*PI*i/k);
                auto cur_err=std::abs(f(c));
                if(cur_err<err)
                {
                    err=cur_err;
                    ret=c;
                }
            }
//            std::cout<<ret<<std::endl;
            return ret;
        }
        std::vector<coef_type> roots()const
        {
            int n=this->degree();
            if(n<=0)
                throw std::runtime_error("degree is less than zero");
            if(n==1)
                return {-(*this)[0]/(*this)[1]};
            if(n==2)
            {
                auto a=(*this)[2],b=(*this)[1],c=(*this)[0];
                auto d=std::sqrt(b*b-base_type(4)*(a*c)),t=(0.5)/a;
                return {(d-b)*t,(-d-b)*t};
            }
            auto f=*this;
            std::vector<coef_type> ans(n);
            for(int i=0,guard=n-2;i<guard;++i)
            {
                auto temp=f.root_without_init();
                ans[i]=temp;
                f=f.div_monomial_factor(temp).first;
            }
            auto a=f[2],b=f[1],c=f[0];
            auto d=std::sqrt(b*b-base_type(4)*(a*c)),t=(0.5)/a;
            ans[n-2]=(d-b)*t;
            ans[n-1]=(-d-b)*t;
            for(auto &x:ans)
                x=this->root_on_guess(x);
            return ans;
        }
        cpolynomial mul_monomial_factor(const coef_type &a)const
        {
            return this->basic_polynomial<coef_type>::mul_monomial_factor(a);
        }
        cpolynomial &operator%=(const cpolynomial &v)
        {
            this->basic_polynomial<coef_type>::operator%=(v);
            return *this;
        }
        cpolynomial &operator*=(const cpolynomial<T> &rhs)
        {
            this->basic_polynomial<coef_type>::operator*=(rhs);
            return *this;
        }
        cpolynomial &operator+=(const cpolynomial<T> &rhs)
        {
            this->basic_polynomial<coef_type>::operator+=(rhs);
            return *this;
        }
        cpolynomial &operator-=(const cpolynomial<T> &rhs)
        {
            this->basic_polynomial<coef_type>::operator-=(rhs);
            return *this;
        }
        cpolynomial &operator/=(const cpolynomial<T> &rhs)
        {
            this->basic_polynomial<coef_type>::operator/=(rhs);
            return *this;
        }
        cpolynomial &operator/=(const coef_type &c)
        {
            this->basic_polynomial<coef_type>::operator/=(c);
            return *this;//类型不完全匹配时，调用了构造函数，返回了右值引用
        }
        cpolynomial &operator*=(const coef_type &c)
        {
            this->basic_polynomial<coef_type>::operator*=(c);
            return *this;
        }
        std::pair<cpolynomial,cpolynomial> polydiv(const cpolynomial<T> &rhs)const
        {
            return this->basic_polynomial<coef_type>::polydiv(rhs);
        }
        cpolynomial polymul(const cpolynomial<T> &rhs)const
        {
            return this->basic_polynomial<coef_type>::polymul(rhs);
        }
        cpolynomial polyadd(const cpolynomial<T> &rhs)const
        {
            return this->basic_polynomial<coef_type>::polyadd(rhs);
        }
        cpolynomial polysub(const cpolynomial<T> &rhs)const
        {
            return this->basic_polynomial<coef_type>::polysub(rhs);
        }
        cpolynomial derivate()const
        {
            return this->basic_polynomial<coef_type>::derivate();
        }
        std::pair<cpolynomial,coef_type> div_monomial_factor(const coef_type &a)const
        {
            return this->basic_polynomial<coef_type>::div_monomial_factor(a);
        }
    };

    template<class T>
    inline cpolynomial<T> operator+(const cpolynomial<T> &lhs,const cpolynomial<T> &rhs)
    {
        return lhs.polyadd(rhs);
    }

    template<class T>
    inline cpolynomial<T> operator-(const cpolynomial<T> &lhs,const cpolynomial<T> &rhs)
    {
        return lhs.polysub(rhs);
    }

    template<class T>
    inline cpolynomial<T> operator/(const cpolynomial<T> &lhs,const cpolynomial<T> &rhs)
    {
        return lhs.polydiv(rhs).first;
    }

    template<class T>
    inline cpolynomial<T> operator%(const cpolynomial<T> &lhs,const cpolynomial<T> &rhs)
    {
        return lhs.polydiv(rhs).second;
    }

    template<class T>
    inline cpolynomial<T> operator*(const cpolynomial<T> &lhs,const cpolynomial<T> &rhs)
    {
        return lhs.polymul(rhs);
    }
}

#endif // POLYNOMIAL_H
