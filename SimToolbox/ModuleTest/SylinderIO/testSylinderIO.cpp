#include "../../Trilinos/TpetraUtil.hpp"

#include "../../Sylinder/Sylinder.hpp"
#include "../../Util/EigenDef.hpp"

void writeVTK(const std::vector<Sylinder> &sylinder, const std::string &baseFolder) {
    int snapID = 0;
    auto commRcp = getMPIWORLDTCOMM();

    Sylinder::writeVTP(sylinder, baseFolder, std::to_string(snapID), commRcp->getRank());
    // write parallel head
    if (commRcp->getSize() > 1 && commRcp->getRank() == 0) {
        Sylinder::writePVTP(baseFolder, "000", commRcp->getSize());
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    std::vector<Sylinder> sylinder;

    for (int i = 0; i < 10; i++) {
        sylinder.emplace_back(i, 1.0, 1.0, 5.0, 5.0, 50 * Evec3::Random(), Equatn::UnitRandom());
    }

    writeVTK(sylinder, "./");

    MPI_Finalize();

    return 0;
}