#ifndef NBNODE_HPP
#define NBNODE_HPP

#include <mpi.h>
#include <omp.h>

#include <algorithm>
#include <cstdio>
#include <deque>
#include <iostream>
#include <tuple>
#include <type_traits>
#include <vector>

#include "Buffer.hpp"
#include "MortonID.hpp"

#define MAX_DEPTH_NB 15
// current  implementation 3D only
#define COLLEAGUECOUNT 27

template <typename ObjT1, typename ObjT2> // e.g. MyClass<T1, T2>
class NBNode {
    // the node for octree, holding multiple types of objects

  protected:
    size_t node_id;     // For translating node pointer to index.
    int dim = 3;        // Dimension of the tree, should be always 3
    const int maxDepth; // max tree depth
    int status;
    bool ghost;       // ghost node owned by other process
    int depth;        // Depth of the node (root -> 0)
    int path2node;    // Identity among siblings
    NBNode *parent;   // Pointer to parent node
    NBNode *child[8]; // Pointer child nodes, maximum 8 pointers

    size_t maxPts;    // max number of points in all its descendent
    long long weight; // weight for load balancing across mpi ranks
    // weight should be the same as used in PVFMM to generate octree with similar structure

    double coord[3]; // the coord (x,y,z) of this node TODO: check coord is center or the vertex with lower xyz
    NBNode<double> *colleague[COLLEAGUECOUNT]; // pointer to its neighbor nodes.

    // length scale and coordinate base should be the same for all nodes in a tree across all mpi ranks
    // double lengthScale; // the octree is [0,1)^3, multiply this to get the real length
    // double coordBase[3]; // the octree is [0,1)^3, add this to get the real coordinate

    std::vector<ObjT1> *objT1ContainerPtr; // pointer to container of objects
    std::vector<size_t> objT1Index;        // indices of objs owned by this node (non-empty only if it is a leaf node)
    std::vector<double> objT1Coord;        // the coordinate of objs in the first container, used for tree subdivision
    std::vector<size_t> objT1scatter;      // scatter index

    std::vector<ObjT2> *objT2ContainerPtr =
        nullptr;                      // pointer to container of objects of Type 2, = nullptr if T2 is void
    std::vector<size_t> objT2Index;   // empty if T2 is void
    std::vector<size_t> objT2scatter; // scatter index

    Buffer packedData; // packed data for transfer over MPI

  public:
    // constructor
    NBNode(int maxDepth_, int maxPts_, std::vector<ObjT1> *objT1ContainerPtr_,
           std::vector<ObjT2> *objT2ContainerPtr_ = nullptr)
        : maxDepth(maxDepth_ < MAX_DEPTH_NB ? maxDepth_ : MAX_DEPTH_NB), maxPts(maxPts_), parent(nullptr),
          child(nullptr), status(1), ghost(false), depth(0), path2node(0), weight(0), coord({0, 0, 0}),
          colleague(nullptr), node_id(0) {

        objT1ContainerPtr = objT1ContainerPtr_;
        if (std::is_void<ObjT2>::value == true) {
            objT2ContainerPtr = nullptr;
        } else {
            // objT2ContainerPtr = objT2ContainerPtr_;
            // TODO: octree with two species is not implemented yet.
            printf("Octree with two species is not implemented yet.");
            exit(1);
        }
    }

    // destructor. TODO: destruct an octant without invalidating its parent's and children's pointers to it
    ~NBNode() {
        // delete all children
        for (int i = 0; i < 8; i++) {
            if (child[i] != nullptr) {
                delete child[i];
                child[i] = nullptr;
            }
        }
    }

    // forbit copy
    NBNode(const NBNode &) = delete;
    NBNode(NBNode &&) = delete;
    const NBNode &operator=(const NBNode &) = delete;
    const NBNode &operator=(NBNode &&) = delete;

    // Initialize the node by passing the relevant data.
    void initialize(NBNode *parent_, int path2node_) {
        // TODO:  check this
        parent = parent_;
        depth = (parent == nullptr ? 0 : parent->getDepth() + 1);
        if (parent != nullptr) {
            dim = parent->dim;
            maxDepth = parent->maxDepth;
            maxPts = parent->maxPts;
        }
        setPath2Node(path2node_);
    }

    void clearData() {
        // TODO:
    }

    // Returns the dimension of the tree.
    int getDim() { return dim; }

    // Returns the depth of this node. (Root has depth 0)
    int getDepth() { return depth; }

    // Returns 'true' if this is a leaf node.
    bool isLeaf() { return (child == nullptr); }

    // Returns the child corresponding to the input parameter.
    NBNode *getChild(int id) {
        assert(id < (1 << dim));
        if (child == nullptr) {
            return nullptr;
        } else {
            return child[id];
        }
    }

    // Returns a pointer to the parent node.
    NBNode *getParent() { return parent; };

    /**
     * Returns the index which corresponds to this node among its
     * siblings (parent's children).
     * this->getParent()->getChild(this->getPath2Node())==this
     */
    int getPath2Node() { return path2node; }

    void setPath2Node(int path2node_) {
        assert(path2node_ > 0 && path2node_ < static_cast<int>(1U << dim));
        path2node = path2node_;
    }

    /**
     * Allocate a new object of the same type (as the derived class) and
     * return a pointer to it type cast as (TreeNode*).
     */
    NBNode *allocNewNode(NBNode *n_ = nullptr) {
        NBNode *n =
            (n_ == nullptr ? new NBNode(this->maxDepth, this->maxPts, this->objT1ContainerPtr, this->objT2ContainerPtr)
                           : n_);
        return n;
    }

    /**
     * \brief Evaluates and returns the subdivision condition for this node.
     * 'true' if node requires further subdivision.
     */
    bool getSubdivCond() {
        if (!isLeaf()) {
            int n = (1UL << dim);
            for (int i = 0; i < n; i++) {
                TreeNode *ch = this->Child(i);
                assert(ch != NULL); // This should never happen
                if (!ch->IsLeaf())
                    return true;
            }
            if (Depth() >= max_depth)
                return false;
            return true;
        } else {
            if (this->Depth() >= max_depth)
                return false;
            return false;
        }
    }

    // Create child nodes and Initialize them.
    void subdivide() {
        // TODO: alloc new child nodes, pass data and pointers, set depth, coord, and morton id.
    }

    // Truncates the tree i.e. makes this a leaf node.
    void truncate();

    // Set the parent for this node.
    void setParent(NBNode *p, int path2node_);

    // Set a child for this node.
    void setChild(NBNode *c, int id);

    // Returns status.
    int &getStatus();

    // Update status for all nodes up to the root node.
    void setStatus(int flag);

    /**
     * Returns list of coordinate and value vectors which need to be
     * sorted and partitioned across MPI processes and the scatter index is
     * saved.
     */
    void getNodeDataVec(std::vector<std::vector<double> *> &coord, std::vector<std::vector<size_t> *> &scatter) {
        coord.push_back(&objT1Coord);
        scatter.push_back(&objT1scatter);
    }

    /**
     * \brief Returns the colleague corresponding to the input index.
     */
    NBNode *getColleague(int index) { return colleague[index]; }

    /**
     * \brief Set the colleague corresponding to the input index.
     */
    void setColleague(NBNode *node_, int index) { colleague[index] = node_; }

    /**
     * \brief Returns the cost of this node. Used for load balancing.
     */
    long long &getNodeCost() { return weight; }

    /**
     * \brief Returns an array of size dim containing the coordinates of the
     * node.
     */
    double *getCoord() {
        assert(coord != nullptr);
        return coord;
    }

    /**
     * \brief Determines if the node is a Ghost node or not.
     */
    bool isGhost() { return ghost; }

    /**
     * \brief Sets the ghost flag of this node.
     */
    void setGhost(bool x) { ghost = x; }

    /**
     * \brief Gets Morton Id of this node.
     */
    inline MortonId getMortonId();

    /**
     * \brief Sets the coordinates of this node using the given Morton Id.
     */
    inline void setCoord(MortonId &mid);

    /**
     * \brief Pack this node to be transmitted to another process. The node
     * is responsible for allocating and freeing the memory for the actual data.
     */
    packedData pack(bool ghost = false, size_t offset = 0);

    /**
     * \brief Initialize the node with data from another process.
     */
    void unpack(PackedData data, bool own_data = true);

    /**
     * \brief Read source distribution at points on a grid defined by array of x,
     * y and z coordinates.
     */
    void readVal(std::vector<Real_t> x, std::vector<Real_t> y, std::vector<Real_t> z, Real_t *val,
                 bool show_ghost = true);

    /**
     * \brief Append node VTU data to vectors.
     */
    template <class VTUData_t, class Node_t>
    static void VTU_Data(VTUData_t &vtu_data, std::vector<Node_t *> &nodes, int lod);
};

#endif