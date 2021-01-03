#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include "numeric_type_helper.h"
#include <vector>
#include <cmath>

namespace song
{
    template<class T>
        requires imp::numeric_type_helper<T>::floating_point||
    imp::numeric_type_helper<T>::complex_floating_point
    class matrix
    {
    public:
        typedef typename imp::numeric_type_helper<T>::argument_type element_type;
        typedef typename imp::numeric_type_helper<T>::real_type base_type;
        typedef typename imp::numeric_type_helper<T>::complex_type eigenvalue_type;
        static constexpr auto floating_point=imp::numeric_type_helper<T>::floating_point;
        static constexpr auto complex_floating_point=imp::numeric_type_helper<T>::complex_floating_point;
        static constexpr auto eps=std::sqrt(std::numeric_limits<base_type>::epsilon());//epsilon()==2.22045e-16;
    private:
        std::size_t rows,cols;
        std::vector<element_type> m_data;
    public:
        matrix()=default;
        matrix(int height,int width):rows(height),cols(width),m_data(height*width){}
        matrix(int height,int width,const std::vector<element_type> &data):rows(height),cols(width),m_data(data)
        {
            m_data.resize(rows*cols);
        }
        matrix(int height,int width,std::vector<element_type> &&data):rows(height),cols(width),m_data(std::move(data))
        {
            m_data.resize(rows*cols);
        }
        matrix(const std::vector<element_type> &p):
            std::vector<element_type>(p),rows(this->size()),cols(1){}//实例化要补全类型参数
        matrix(std::vector<element_type> &&p):
            std::vector<element_type>(std::move(p)),rows(this->size()),cols(1){}
        static constexpr auto conj(element_type x)
        {
            if constexpr(floating_point)
                return x;
            else
                return std::conj(x);
        }
        auto &operator()(int i,int j)
        {
            return m_data.at(i*cols+j);
        }
        auto &operator()(int i,int j)const
        {
            return m_data.at(i*cols+j);
        }
        auto resize(int m,int n=1)
        {
            rows=m;
            cols=n;
            return this->std::vector<T>::resize(m*n);
        }
        matrix<element_type> Hessenberg()const
        {
            if(rows!=cols)
                return {};
            int n=rows;
            matrix hessen=*this;
            if(n<=2)
                return hessen;
            for(int l=0;l<n-2;++l)
            {
                std::vector<element_type> u(n-l-1);
                for(int i=l+1,guardi=n;i<guardi;++i)
                    u[i-l-1]=hessen(i,l);
                auto radius=[](const std::vector<element_type> &u)
                {
                    base_type norm=0.0;
                    for(auto &x:u)
                    {
                        if constexpr(floating_point)
                        {
                            norm+=x*x;
                        }
                        else
                        {
                            norm+=std::norm(x);
                        }
                    }
                    return std::sqrt(norm);
                };
                //|x|
                auto norm=radius(u);
                if(norm<eps)
                    continue;
                if constexpr(floating_point)
                {
                    if(u[0]<0)
                        norm=-norm;
                }
                else
                {
                    if(u[0].real()<0)
                        norm=-norm;
                }
                //x+|x|e1
                u[0]+=norm;
                hessen(l+1,l)=-norm;
                for(int i=l+2;i<n;++i)
                    hessen(i,l)=0.0;
                norm=radius(u);
                //u
                for(auto &x:u)
                    x/=norm;
                for(int i=l+1;i<n;++i)
                {
                    //2*uH*api
                    element_type temp=0;
                    for(int j=l+1;j<n;++j)
                    {
                        temp+=conj(u[j-(l+1)])*hessen(j,i);
                    }
                    temp*=base_type(2);
                    //api-u*(2*uH*api)
                    for(int j=l+1;j<n;++j)
                        hessen(j,i)-=u[j-(l+1)]*temp;
                }
                for(int j=l;j<n;++j)
                    hessen(l+1,j)=-hessen(l+1,j);
                 for(int i=l;i<n;++i)
                {
                    //2*apiH*u
                    element_type temp=0;
                    for(int j=l+1;j<n;++j)
                    {
                            temp+=conj(hessen(i,j))*u[j-(l+1)];
                    }
                    temp*=base_type(2);
                    //api-u*(2*uH*api)
                    for(int j=l+1;j<n;++j)
                        hessen(i,j)-=u[j-(l+1)]*temp;
                }
                for(int i=l;i<n;++i)
                    hessen(i,l+1)=-hessen(i,l+1);
            }
            return hessen;
        }
        std::pair<matrix,matrix> HessenQR()const
        {
            if(rows!=cols)
                return {};
            int n=rows;
            if(n<=1)
                return {{1,1,{1.0}},*this};
            if(n<=2)
                return {{2,2,{1.0,0.0,0.0,1.0}},*this};
            matrix Q(n,n),QHA=*this;
            for(int i=0,guardi=n;i<guardi;++i)
            {
                Q(i,i)=element_type(1.0);
            }
            for(int k=0,guardk=n-1;k<guardk;++k)
            {
                auto x1=QHA(k,k),x2=QHA(k+1,k);
                base_type r;
                if constexpr(floating_point)
                {
                    r=std::hypot(x1,x2);
                }
                else
                {
                   r=std::sqrt(std::norm(x1),std::norm(x2));
                }
                if(r<eps)
                    continue;
                auto r11=x1/r,r12=x2/r,r21=-r12,r22=r11;
                //update Q
                matrix Qsub(k+2,2);
                for(int j=0,guardj=k+1;j<guardj;++j)
                {
                    auto q1=Q(j,k);
                    Qsub(j,0)=conj(r11)*q1;
                    Qsub(j,1)=conj(r21)*q1;
                }
                {
                    auto q2=Q(k+1,k+1);
                    Qsub(k+1,0)=conj(r12)*q2;
                    Qsub(k+1,1)=conj(r22)*q2;
                }
                for(int j=0,guardj=k+2;j<guardj;++j)
                {
                    Q(j,k)=Qsub(j,0);
                    Q(j,k+1)=Qsub(j,1);
                }
                //update QTA
                 matrix QHAsub(2,n-k);
                {
                    QHAsub(0,0)=r;
                    QHAsub(1,0)=0;
                }
                for(int j=k+1,guardj=n;j<guardj;++j)
                {
                    auto q1=QHA(k,j),q2=QHA(k+1,j);
                    QHAsub(0,j-k)=r11*q1+r12*q2;
                    QHAsub(1,j-k)=r21*q1+r22*q2;
                }
                for(int j=k,guardj=n;j<guardj;++j)
                {
                    QHA(k,j)=QHAsub(0,j-k);
                    QHA(k+1,j)=QHAsub(1,j-k);
                }
            }
            return {Q,QHA};
        }
        std::vector<element_type> QRIter()const
        {
            if(cols!=rows)
                return {};
            int n=cols;
            auto M=Hessenberg();
            for(int it=0;it<101;++it)
            {
                auto [Q,R]=M.HessenQR();
                for(int i=0;i<n;++i)
                    for(int j=0;j<n;++j)
                    {
                        M(i,j)=element_type(0);
                        for(int k=i,guardk=std::min(j+2,n);k<guardk;++k)
                            M(i,j)+=R(i,k)*Q(k,j);
                    }
                n=cols;
            }
            return {};
        }
    };
}



#endif // MATRIX_H_INCLUDED
