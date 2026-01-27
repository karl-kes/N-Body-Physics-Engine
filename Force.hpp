#include <vector>

class Particle;

class Force_Law {
public:
    virtual ~Force_Law() = default;
    virtual void force( std::vector<Particle> particles ) = 0;
};

class Gravity : public Force_Law {
public:
    Gravity();
    void force( std::vector<Particle> particles ) override {

    } 
};