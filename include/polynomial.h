#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <vector>
#include <cmath>
#include <complex>
#include <algorithm>
#include <type_traits>
#include <concepts>
#include <numbers>

namespace song
{
    namespace imp
    {
        template<class T>
        struct polynomial_type_helper
        {
            typedef T coefficient_type;
            typedef T base_type;
            typedef std::complex<T> root_type;
            static constexpr bool floating_point=std::floating_point<T>;
            static constexpr bool complex_floating_point=false;
        };
        template<class T>
        struct polynomial_type_helper<std::complex<T>>
        {
            typedef std::complex<T> coefficient_type;
            typedef T base_type;
            typedef std::complex<T> root_type;
            static constexpr bool floating_point=false;
            static constexpr bool complex_floating_point=std::floating_point<T>;
        };
    }
    template<class T>
    requires imp::polynomial_type_helper<T>::floating_point||
        imp::polynomial_type_helper<T>::complex_floating_point
    class polynomial:public std::vector<T>
    {
    public:
        typedef typename imp::polynomial_type_helper<T>::coefficient_type coefficient_type;
        typedef typename imp::polynomial_type_helper<T>::base_type base_type;
        typedef typename imp::polynomial_type_helper<T>::root_type root_type;
        static constexpr auto floating_point=imp::polynomial_type_helper<T>::floating_point;
        static constexpr auto complex_floating_point=imp::polynomial_type_helper<T>::complex_floating_point;
        static constexpr auto eps=std::numeric_limits<base_type>::epsilon();//epsilon()==2.22045e-16;
        static constexpr auto inf=std::numeric_limits<base_type>::infinity();
        static constexpr auto PI=std::numbers::pi_v<base_type>;
    public:
        coefficient_type root_with_init(const coefficient_type &arg_x)const
        requires complex_floating_point
        {
            auto x=arg_x;
            auto dx=offset(x);
            auto err=std::abs((*this)(x));
            while(err>eps)
            {
                auto cur_x=x+dx;
                auto err1=std::abs((*this)(cur_x));
                bool finish_iter=false;
                if(err1>err)
                {
                    while(true)
                    {
                        if(std::abs(dx)<eps)
                        {
                            finish_iter=true;
                            break;
                        }
                        auto dx2=dx;
                        auto err2=err1;
                        for(int i=1;i<=2;++i)
                        {
                            auto dx3=dx*std::polar(base_type(1.0),i*(PI*base_type(2.0/3.0)));
                            auto cur_x3=x+dx3;
                            auto err3=std::abs((*this)(cur_x3));
                            if(err3<err2)
                            {
                                dx2=dx3;
                                err2=err3;
                            }
                        }
                        if(err>err2)
                        {
                            cur_x=x+dx2;
                            err1=err2;
                            break;
                        }
                        dx*=base_type(0.5);
                    }
                }
                if(finish_iter)
                    break;
                x=cur_x;
                dx=offset(x);
                err=err1;
            }
            return x;
        }
        root_type root_without_init()const
        {
            int deg=degree();
            if(deg<=0)
                throw std::runtime_error("degree is less than zero");
            auto off=-(*this)[deg-1]/(coefficient_type(deg));
//            auto tempf=translation(off);
            coefficient_type a=(*this)(off);
            base_type h=1.0/deg;
            base_type x_r=std::pow(std::abs(a),h);
            coefficient_type x=off;
            base_type err=std::abs((*this)(x));
            for(int i=0;i<deg;++i)
            {
                auto curx=off+std::polar(x_r,2*i*(PI*h));
                auto curerr=std::abs((*this)(curx));
                if(curerr<err)
                {
                    x=curx;
                    err=curerr;
                }
            }
            return root_with_init(x);
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
                return {inf};
            auto fz=(*this)(z),gz=g(z);
            if(std::abs(gz)>eps)
                return fz/gz;
            if(std::abs(fz)>eps)
                return {inf};
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
                    return {inf};
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
            int nc=degree(),nd=nc;
            polynomial pd(this->size());
            pd[0]=(*this)[nc];
            for(int i=nc-1;i>=0;--i)
            {
                int nnd=std::min(nd,nc-i);
                for(int j=nnd;j>=1;--j)
                    pd[j]=pd[j]*x+pd[j-1];
                pd[0]=pd[0]*x+(*this)[i];
            }
            return pd;
        }

        std::vector<root_type> roots_of_degree1()const
        {
            return {-root_type((*this)[0])};
        }
        std::vector<root_type> roots_of_degree2()const
        {
            auto b=(*this)[1],c=(*this)[0];
            auto delta=b*b-base_type(4.0)*c;
            if constexpr(floating_point)
            {
                if(delta>=0)
                {
                     coefficient_type sqrt_delta=std::sqrt(delta);
                    auto q=-0.5*(b+std::copysign(sqrt_delta,b));
                    return {q,c/q};
                }
                else
                {
                    auto sqrt_delta=std::sqrt(-delta);
                    return {root_type(-0.5*b,-0.5*sqrt_delta),root_type(-0.5*b,0.5*sqrt_delta)};
                }
            }
            else
            {
                auto sqrt_delta=std::sqrt(b*b-base_type(4.0)*c);
                auto Re=(b*sqrt_delta).real();
                if(Re<=0)
                    sqrt_delta=-sqrt_delta;
                auto qq=-base_type(0.5)*(b+sqrt_delta);
                return {qq,c/qq};
            }
            return {};
        }
        //from https://zhuanlan.zhihu.com/p/40349993
        std::vector<root_type> roots_of_degree3()const
        {
            auto a=(*this)[2],b=(*this)[1],c=(*this)[0];
            auto a2=a*a,a3=a2*a;
            root_type Q=(a2-base_type(3.0)*b)/base_type(9.0),R=(base_type(2.0)*a3-base_type(9.0)*a*b+base_type(27.0)*c)/base_type(54.0);
            auto Q1d2=std::sqrt(Q),Q3d2=Q*Q1d2;
            if(std::abs(R)<std::abs(Q3d2))
            {
                auto theta=std::acos(R/Q3d2);
                auto x1=-base_type(2.0)*Q1d2*std::cos(theta/base_type(3.0))-a/base_type(3.0);
                auto x2=-base_type(2.0)*Q1d2*std::cos((theta+2*PI)/base_type(3.0))-a/base_type(3.0);
                auto x3=-base_type(2.0)*Q1d2*std::cos((theta-2*PI)/base_type(3.0))-a/base_type(3.0);
                return {x1,x2,x3};
            }
            auto R2=R*R,Q3=Q3d2*Q3d2,Temp=std::sqrt(R2-Q3);
            if(std::real(R*Temp)<0)
                Temp=-Temp;
            auto A=std::pow(R+Temp,base_type(1.0/3));
            auto B=std::abs(A)>eps?Q/A:base_type(0.0);
            auto x1=(A+B)-a*base_type(1.0/3.0);
            auto x2=-base_type(0.5)*(A+B)-a*base_type(1.0/3.0)+std::sin(PI/3)*(A-B)*root_type(0.0,1.0);
            auto x3=-base_type(0.5)*(A+B)-a*base_type(1.0/3.0)-std::sin(PI/3)*(A-B)*root_type(0.0,1.0);
            return {x1,x2,x3};
        }
        //from https://www.cnblogs.com/larissa-0464/p/11706131.html
        std::vector<root_type> roots_of_degree4()const
        {
            auto b=(*this)[3];
            auto c=(*this)[2];
            auto d=(*this)[1];
            auto e=(*this)[0];
            auto c2=c*c,bd=b*d,b2=b*b;
            auto P=(c2+base_type(12.0)*e-base_type(3.0)*bd)/base_type(9.0);
            auto Q = (base_type(27.0)*d*d+base_type(2.0)*c2*c+base_type(27.)*b2*e-base_type(72.)*c*e-base_type(9.)*bd*c)/base_type(54.);
            root_type D = std::sqrt(Q*Q-P*P*P);
            if(std::real(Q*D)<0)
                D=-D;
            auto t1=Q+D;
            auto u=std::pow(t1,base_type(1.0/3));
            auto v=std::abs(u)<eps?base_type(0.0):P/u;
            root_type w[]={std::polar<base_type>(1.,2./3*PI),std::polar<base_type>(1.,-2./3*PI)};
            auto temp0=b2-base_type(8./3)*c;
            auto temp1=base_type(4.)*(u+v);
            auto sqr_m=temp0+temp1;
            auto abs_sqr_m=std::abs(sqr_m);
            for(int i=0;i<2;++i)
            {
                auto tmp1=base_type(4.)*(w[i]*u+w[1-i]*v);
                auto sqr_temp=temp0+tmp1;
                auto abs_sqr_temp=std::abs(sqr_temp);
                if(abs_sqr_m<abs_sqr_temp)
                {
                    sqr_m=sqr_temp;
                    abs_sqr_m=abs_sqr_temp;
                    temp1=tmp1;
                }
            }
            auto m=std::sqrt(sqr_m);
            coefficient_type S=temp0,TT=0.0;
            if(std::abs(m)>eps)
            {
                S=base_type(2.)*temp0-temp1;
                TT=(base_type(8.)*b*c-base_type(16.)*d-base_type(2.)*b2*b)/m;
            }
            auto bm0=-b-m,bm1=-b+m;
            auto sq_ST0=std::sqrt(S-TT),sq_ST1=std::sqrt(S+TT);
            return {(bm0+sq_ST0)/base_type(4.),(bm0-sq_ST0)/base_type(4.),(bm1+sq_ST1)/base_type(4.),(bm1-sq_ST1)/base_type(4.)};
        }
        coefficient_type offset(const coefficient_type &x)const
        {
            auto dpn=derivate(),d2pn=dpn.derivate();
            auto pnx=(*this)(x),dpnx=dpn(x),d2pnx=d2pn(x);
            int n=degree();
            auto flag1=std::abs(dpnx)<eps,flag0=std::abs(pnx)<eps;
            if(!flag1)
            {
                auto inv_dpnx=base_type(1.0)/dpnx,temp1=pnx*inv_dpnx,temp2=temp1*d2pnx*inv_dpnx;
                auto delta1d2=std::sqrt(base_type(n-1)*(base_type(n-1)-base_type(n)*temp2));
                if(delta1d2.real()<0)
                    delta1d2=-delta1d2;
                auto d1=base_type(1.0)+delta1d2;
                return -temp1*(base_type(n)/d1);
            }
            // flag1
            if(!flag0)
            {
                //dpnx==0,pnx!=0,
                if(std::abs(d2pnx)<eps)
                    return {inf,inf};
                auto delta=std::sqrt(base_type(n-1)*(base_type(n-1)*(dpnx*dpnx)-base_type(n)*(pnx*d2pnx)));
                if((dpnx*delta).real()<0)
                    delta=-delta;
                auto d1=dpnx+delta;
                return -base_type(n)*(pnx/d1);
            }
            // flag1&&flag2
            auto temp1=div_on_point(dpn,x);
            auto temp2=((*this)*d2pn).div_on_point(dpn*dpn,x);
            return -temp1/(base_type(1.0)-temp2);
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
        requires complex_floating_point
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

        polynomial &operator%=(const polynomial &v)
        {
            *this=poly_divide(v).second;
            return *this;
        }
        polynomial &operator*=(const polynomial &rhs);

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
            return *this;
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
            return *this;
        }
        polynomial &operator/=(const polynomial &v)
        {
            *this=poly_divide(v).first;
            return *this;
        }
        polynomial &operator*=(const coefficient_type &c)
        {
            auto tempc=c;
            for(auto &pi:*this)
                pi*=tempc;
            return *this;
        }
        polynomial &operator/=(const coefficient_type &c)
        {
            return (*this)*=base_type(1.0)/c;
        }

        //规范化，不允许0系数占据最高次
        polynomial normalize()const
        {
            auto p=*this;
            for(int i=degree();i>=0;--i)
                if(std::abs(p[i])>eps)
                {
                   p.erase(this->cbegin()+i+1,this->cend());
                    break;
                }
            return p;
        }
        polynomial monic()const
        {
            polynomial ret(this->size());
            auto inv_an=coefficient_type(1.0)/this->back();
            for(int i=0,guard=this->size();i<guard;++i)
            {
                ret[i]=(*this)[i]*inv_an;
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
            return std::pair{std::move(q),std::move(r)};
        }

        polynomial derivate()const
        {
            int n=this->size();
            if(n<=1)
                return {};
            polynomial p(n-1);
            for(int i=1;i<n;++i)
            {
                p[i-1]=base_type(i)*(*this)[i];
            }
            return p;
        }
    };

    template<class T>
    inline polynomial<T> operator+(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        auto m=lhs.size(),n=rhs.size();
        if(m<n)
            return rhs+lhs;
        auto p=lhs;
        p+=rhs;
        return p;
    }

    template<class T>
    inline polynomial<T> operator+(const polynomial<T> &p)
    {
        return p;
    }

    template<class T>
    inline auto operator-(const polynomial<T> &p)
    {
        auto ret=p;
        for(auto &c:ret)
            c=-c;
        return ret;
    }

    template<class T>
    inline polynomial<T> operator-(const polynomial<T> &lhs,const polynomial<T> &rhs)
    {
        return lhs+(-rhs);
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
        int m=lhs.size(),n=rhs.size();
        int p=m+n-1;
        polynomial<T> c(p);
        for(int i=0;i<m;++i)
            for(int j=0;j<n;++j)
                c[i+j]+=lhs[i]*rhs[j];
        return c;
    }

    template<class T>
    inline polynomial<T> &polynomial<T>::operator*=(const polynomial &rhs)
    {
        *this=(*this)*rhs;
        return *this;
    }

    template<class T>
    inline polynomial<T> from_roots(const std::vector<T> &roots)
    {
        polynomial<T> ret={1.0};
        for(auto &c:roots)
        {
            ret*={-c,1.0};
        }
        return ret;
    }

    template<class T>
    using cpolynomial=polynomial<std::complex<T>>;
}

#endif // POLYNOMIAL_H
