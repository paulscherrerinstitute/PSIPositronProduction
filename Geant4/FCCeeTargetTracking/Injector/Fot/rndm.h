#ifndef RNDM_SEEN
#define RNDM_SEEN

#include <random>

class RNDM {

  std::default_random_engine generator;
  std::uniform_real_distribution<double> uniform_rng;
  std::normal_distribution<double> normal_rng;

 public:
  
 RNDM(int seed = 1 ) : generator(seed) {}

  inline void set_seed (int seed){
    generator.seed(seed);
  }
  
  inline double rndm() 
  {
    return uniform_rng(generator);
  }

  inline double gaussian(double mean=0.0, double stddev=1.0 )
  {
    return mean + normal_rng(generator) * stddev;
  }

};

#endif
