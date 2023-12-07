//~/homebrew/bin/g++-13 -I /Users/mgrecu/frugally-deep/include/ -I /Users/mgrecu/homebrew/include/ main.cpp


#include <fdeep/fdeep.hpp>
#include<stdio.h>

const auto model = fdeep::load_model("fdeep_model.json");

extern "C" void mp_model_(float *v1, float *v2, float *v3, float *v4, float *o1, float *o2)
{
    const auto result = model.predict(
                        {fdeep::tensor(fdeep::tensor_shape(static_cast<std::size_t>(4)),
                        std::vector<float>{*v1, *v2, *v3, *v4})});
        const std::vector<float> vec = result[0].to_vector();
        *o1=vec[0];
        *o2=vec[1];
}
