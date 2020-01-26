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
    public:
        typedef T coefficient_type;
        using abs_type=decltype(std::abs(coefficient_type()));
        static constexpr auto eps=1000*std::numeric_limits<abs_type>::epsilon();//2.22045e-13;
        static constexpr auto inf=std::numeric_limits<abs_type>::infinity();
        static constexpr auto PI=4*std::atan(abs_type(1.0));
    protected:
        using std::vector<coefficient_type>::vector;//继承vector构造函数，不包括从vector到poly的转换
//        template<template<class> class Container>
        basic_polynomial &correct()
        {
            auto i=this->degree();
            for(;i>=0;--i)
            {
                auto x=std::pow(std::abs((*this)[i]),1.0/i);
                if(x>eps)
                    break;
                --i;
            }
            this->erase(this->begin()+i+1,this->end());
            return *this;
        }
        basic_polynomial translation(const coefficient_type &x)const
        {
            auto dfi_div_ifact=*this;
            basic_polynomial ret(this->size());
            ret[0]=dfi_div_ifact(x);
            for(int i=1,guard=ret.size();i<guard;++i)
            {
                dfi_div_ifact=dfi_div_ifact.derivate();
                dfi_div_ifact/=coefficient_type(i);
                ret[i]=dfi_div_ifact(x);
            }
            return ret;
        }
    public:
        basic_polynomial(const std::vector<coefficient_type> &v):std::vector<coefficient_type>((v)){}
        basic_polynomial(std::vector<coefficient_type> &&v):std::vector<coefficient_type>((std::move(v))){}
//        basic_polynomial(std::initializer_list<coefficient_type> il):std::vector<coefficient_type>(fix_container(il)){}
        coefficient_type operator()(const coefficient_type &x)const
        {
            if(this->empty())
                return coefficient_type(0.0);
            const int n=this->size()-1;
            coefficient_type p=(*this)[n];
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
        basic_polynomial mul_monomial_factor(const coefficient_type &a)const
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
        auto div_monomial_factor(const coefficient_type &a)const
        //f/(x-a)
        {
            int n=this->size()-1;
            coefficient_type rem=(*this)[n];
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
            int n=this->size()-1,nv=v.size()-1;
            basic_polynomial r=*this,q(n);
            for(int k=n-nv;k>=0;--k)
            {
                q[k]=r[nv+k]/v[nv];
                for(int j=nv+k-1;j>=k;--j)
                    r[j]-=q[k]*v[j-k];
            }
            this->swap(q);
            return *this;
        }
        basic_polynomial &operator%=(const basic_polynomial &v)
        {
            int n=this->size()-1,nv=v.size()-1;
            basic_polynomial q(n);
            for(int k=n-nv;k>=0;--k)
            {
                q[k]=(*this)[nv+k]/v[nv];
                for(int j=nv+k-1;j>=k;--j)
                    (*this)[j]-=q[k]*v[j-k];
            }
            this->erase(this->cbegin()+nv,this->cend());
            this->correct();
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
                p[i-1]=abs_type(i)*(*this)[i];
            }
            return p;
        }
        basic_polynomial &operator*=(const basic_polynomial &rhs)
        {
            return *this=this->polymul(rhs);
        }
        basic_polynomial &operator*=(const coefficient_type &c)
        {
            for(auto &pi:*this)
                pi*=c;
            return *this;
        }
        basic_polynomial &operator/=(const coefficient_type &c)
        {
            return (*this)*=1.0/c;
        }
        int degree()const
        //最高次数，值为-1表示零多项式
        {
            return int(this->size())-1;
        }
    };


    template<class T>
    class polynomial:public basic_polynomial<T>
    {
    public:
        typedef typename basic_polynomial<T>::coefficient_type coefficient_type;
        typedef typename basic_polynomial<coefficient_type>::abs_type abs_type;

        using basic_polynomial<T>::basic_polynomial;
        polynomial(const basic_polynomial<T> &rhs):basic_polynomial<T>(rhs){}
        polynomial(basic_polynomial<T> &&rhs):basic_polynomial<T>(std::move(rhs)){}
    };

    template<class T>
    class polynomial<std::complex<T>>:public basic_polynomial<std::complex<T>>
    {
    public:
        typedef typename basic_polynomial<std::complex<T>>::coefficient_type coefficient_type;
        typedef typename basic_polynomial<coefficient_type>::abs_type abs_type;

//    private:

        coefficient_type root_on_guess(const coefficient_type &arg_x)const
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
        coefficient_type root_without_init()const
        {
            int deg=this->degree();
            if(deg<=0)
                throw std::runtime_error("degree is less than zero");
            coefficient_type a=this->front()/this->back();
            abs_type h=1.0/deg;
            abs_type x_r=std::pow(std::abs(a),h);
            coefficient_type x=0,dx=this->offset(x);
            abs_type err=std::abs((*this)(dx));
            for(int i=0;i<deg;++i)
            {
                auto curx=std::polar(x_r,2*i*(this->PI*h));
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
        coefficient_type div_on_point(const polynomial &g,const coefficient_type &z)const
        {
            if(g.empty())
                throw std::runtime_error("g is zero in polynomial<...>::div_on_point");
            auto fz=(*this)(z),gz=g(z);
            if(std::abs(gz)>this->eps)
                return fz/gz;
            if(std::abs(fz)>this->eps)
                return {this->inf,this->inf};
            auto [q,r]=this->polydiv(g);
            auto dr=r.derivate(),dg=g.derivate();
            while(true)
            {
                auto dgz=dg(z),drz=dr(z);
                if(std::abs(dgz)>this->eps)
                {
                    return q(z)+drz/dgz;
                }
                if(std::abs(drz)>this->eps)
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
        polynomial translation(const coefficient_type &x)const
        {
            return this->basic_polynomial<coefficient_type>::translation(x);
        }
    public:
        using basic_polynomial<coefficient_type>::basic_polynomial;
        polynomial(const basic_polynomial<coefficient_type> &p):basic_polynomial<coefficient_type>(p){}//实例化要补全类型参数
        polynomial(basic_polynomial<coefficient_type> &&p):basic_polynomial<coefficient_type>(std::move(p)){}
        coefficient_type offset(const coefficient_type &x)const
        {
            auto dpn=this->derivate(),d2pn=dpn.derivate();
            auto pnx=(*this)(x),dpnx=dpn(x),d2pnx=d2pn(x);
            int n=this->degree();
            auto flag1=std::abs(dpnx)<this->eps,flag0=std::abs(pnx)<this->eps;
            if(!flag1)
            {
                auto t=1.0/dpnx,temp1=pnx*t,temp2=temp1*d2pnx*t;
                auto delta=std::sqrt(abs_type(n-1)*(abs_type(n-1)-abs_type(n)*temp2));
                auto d1=1.0+delta,d2=1.0-delta;
                return temp1*(abs_type(n)/(std::abs(d1)>std::abs(d2)?d1:d2));
            }
            // flag1
            if(!flag0)
            {
                if(std::abs(d2pnx)<this->eps)
                    return {this->inf,this->inf};
                auto delta=std::sqrt(abs_type(n-1)*(abs_type(n-1)*(dpnx*dpnx)-abs_type(n)*(pnx*d2pnx)));
                auto d1=dpnx+delta,d2=dpnx-delta;
                return abs_type(n)*(pnx/(std::abs(d1)>std::abs(d2)?d1:d2));
            }
            // flag1&&flag2
            auto temp1=this->div_on_point(dpn,x);
            auto temp2=this->polymul(d2pn).div_on_point(dpn.polymul(dpn),x);
            return temp1/(1.0-temp2);
        }
        coefficient_type tiny_offset(const coefficient_type &x)const
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
            coefficient_type ret;
            for(int i=0;i<k;++i)
            {
                auto c=std::polar(d,alpha+2*this->PI*i/k);
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
        std::vector<coefficient_type> roots()const
        {
            int n=this->degree();
            if(n<=0)
                throw std::runtime_error("degree is less than zero");
            if(n==1)
                return {-(*this)[0]/(*this)[1]};
            if(n==2)
            {
                auto a=(*this)[2],b=(*this)[1],c=(*this)[0];
                auto d=std::sqrt(b*b-abs_type(4)*(a*c)),t=(0.5)/a;
                return {(d-b)*t,(-d-b)*t};
            }
            auto f=*this;
            std::vector<coefficient_type> ans(n);
            for(int i=0,guard=n-2;i<guard;++i)
            {
                auto temp=f.root_without_init();
                ans[i]=temp;
                f=f.div_monomial_factor(temp).first;
            }
            auto a=f[2],b=f[1],c=f[0];
            auto d=std::sqrt(b*b-abs_type(4)*(a*c)),t=(0.5)/a;
            ans[n-2]=(d-b)*t;
            ans[n-1]=(-d-b)*t;
            for(auto &x:ans)
                x=this->root_on_guess(x);
            return ans;
        }
        polynomial mul_monomial_factor(const coefficient_type &a)const
        {
            return this->basic_polynomial<coefficient_type>::mul_monomial_factor(a);
        }
        polynomial &operator%=(const polynomial &v)
        {
            this->basic_polynomial<coefficient_type>::operator%=(v);
            return *this;
        }
        polynomial &operator*=(const polynomial &rhs)
        {
            this->basic_polynomial<coefficient_type>::operator*=(rhs);
            return *this;
        }
        polynomial &operator+=(const polynomial &rhs)
        {
            this->basic_polynomial<coefficient_type>::operator+=(rhs);
            return *this;
        }
        polynomial &operator-=(const polynomial &rhs)
        {
            this->basic_polynomial<coefficient_type>::operator-=(rhs);
            return *this;
        }
        polynomial &operator/=(const polynomial &rhs)
        {
            this->basic_polynomial<coefficient_type>::operator/=(rhs);
            return *this;
        }
        polynomial &operator/=(const coefficient_type &c)
        {
            this->basic_polynomial<coefficient_type>::operator/=(c);
            return *this;//类型不完全匹配时，调用了构造函数，返回了右值引用
        }
        polynomial &operator*=(const coefficient_type &c)
        {
            this->basic_polynomial<coefficient_type>::operator*=(c);
            return *this;
        }
        std::pair<polynomial,polynomial> polydiv(const polynomial &rhs)const
        {
            return this->basic_polynomial<coefficient_type>::polydiv(rhs);
        }
        polynomial polymul(const polynomial &rhs)const
        {
            return this->basic_polynomial<coefficient_type>::polymul(rhs);
        }
        polynomial polyadd(const polynomial &rhs)const
        {
            return this->basic_polynomial<coefficient_type>::polyadd(rhs);
        }
        polynomial polysub(const polynomial &rhs)const
        {
            return this->basic_polynomial<coefficient_type>::polysub(rhs);
        }
        polynomial derivate()const
        {
            return this->basic_polynomial<coefficient_type>::derivate();
        }
        std::pair<polynomial,coefficient_type> div_monomial_factor(const coefficient_type &a)const
        {
            return this->basic_polynomial<coefficient_type>::div_monomial_factor(a);
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
        return lhs.polydiv(rhs).first;
    }

    template<class T>
    inline polynomial<T> operator%(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.polydiv(rhs).second;
    }

    template<class T>
    inline polynomial<T> operator*(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs.polymul(rhs);
    }

    template<class T>
    basic_polynomial<T> from_roots(const std::vector<T> &roots)
    {
        basic_polynomial<T> ret={1.0};
        for(auto &c:roots)
        {
            ret=ret.mul_monomial_factor(c);
        }
        return ret;
    }
}

#endif // POLYNOMIAL_H
