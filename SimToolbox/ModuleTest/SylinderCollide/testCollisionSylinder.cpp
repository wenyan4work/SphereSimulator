#include "../../Collision/CollisionSylinder.hpp"

int main() {
    // Sylinder sy1(1, 1.0, 1.5, 5.0, 5.0, 3 * Evec3::Random(), Equatn::UnitRandom());
    // Sylinder sy2(2, 1.0, 1.5, 5.0, 5.0, 3 * Evec3::Random(), Equatn::UnitRandom());
    Sylinder sy1(1, 1.0, 1.5, 5.0, 5.0, 3 * Evec3::Random(), Equatn::Identity());
    Sylinder sy2(2, 1.0, 1.5, 5.0, 5.0, 3 * Evec3::Random(), Equatn::Identity());

    CollisionSylinder sycol1, sycol2;
    sycol1.CopyFromFull(sy1);
    sycol2.CopyFromFull(sy2);

    CollisionBlock cb;
    bool col = sycol1.collide(sycol2, cb);
    std::cout << "locI " << cb.posI.transpose() << ", locJ " << cb.posJ.transpose() << std::endl;
    std::cout << "posI " << (sy1.pos).transpose() << ", posJ " << (sy2.pos).transpose()
              << std::endl;
    std::cout << "phi0: " << cb.phi0 << std::endl;

    return 0;
}