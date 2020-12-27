#ifndef NUMERIC_TYPE_HELPER_H_INCLUDED
#define NUMERIC_TYPE_HELPER_H_INCLUDED

#include <complex>

namespace song
{
     namespace imp
    {
        template<class T>
        struct numeric_type_helper
        {
            typedef T argument_type;
            typedef T real_type;
            typedef std::complex<T> complex_type;
            static constexpr bool floating_point=std::floating_point<T>;
            static constexpr bool complex_floating_point=false;
        };
        template<class T>
        struct numeric_type_helper<std::complex<T>>
        {
            typedef std::complex<T> argument_type;
            typedef T real_type;
            typedef std::complex<T> complex_type;
            static constexpr bool floating_point=false;
            static constexpr bool complex_floating_point=std::floating_point<T>;
        };
    }
}

#endif // NUMERIC_TYPE_HELPER_H_INCLUDED
