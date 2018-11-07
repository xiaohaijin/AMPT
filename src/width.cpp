#include <cmath>

using namespace std;

extern "C"{
  float width_(float const &dmass);
}

float pow2(float value);
float pow3(float value);

float width_(float const& dmass)
{
    float return_value = 0.f;

    
    float avmass = 0.938868f;
    float pimass = 0.137265f;
    float qavail = 0.f;

    float aux = 0.25f * pow2((pow2(dmass) - pow2(avmass) - pow2(pimass)))
                - pow2((avmass * pimass));
    if (aux > 0.f) {
        qavail = sqrt(aux / pow2(dmass));
    } else {
        qavail = 1.e-06f;
    }
    return_value = 0.47f * pow3(qavail)
                   / (pow2(pimass) * (1.f + 0.6f * pow2((qavail / pimass))));

    return return_value;
}

float pow2(float value)
{
    return value * value;
}

float pow3(float value)
{
    return value * value * value;
}

