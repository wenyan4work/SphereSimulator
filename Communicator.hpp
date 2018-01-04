#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

#include <deque>
#include <vector>

#include "Buffer.hpp"
#include "ZDD.hpp"

template <class ContainerIn, class ContainerOut, class Object>
class Communicator {
    // construct with a given communication pattern, and then reuse the pattern
    // input to constructor: a const pointer to a const container of objects. a list of destination ranks.
    // input to constructor: a const pointer to a container to hold the received objects

    // objects should define pack() and unpack() to convert between raw char arrays and objects
    // each object should have the same type but variable length
    // each object should have a positive and unique integer as their global id (gid)

  public:
    // constructor
    Communicator(const ContainerIn *const inPtr_, ContainerOut *const outPtr_,
                 const std::vector<std::vector<int>> &sendRanks) noexcept
        : inPtr(inPtr_), outPtr(outPtr_) {}

    // destructor
    ~Communicator() {}

    // forbid copy
    Communicator(const Communicator &other) = delete;
    Communicator(Communicator &&other) = delete;

    Communicator &operator=(const Communicator &) = delete;
    Communicator &operator=(Communicator &&) = delete;

    void transfer(); // the objects are copied and sent to their destination ranks

  private:
    // buffer

    // pointer
    const ContainerIn *const inPtr;
    ContainerOut *const outPtr;
};

#endif