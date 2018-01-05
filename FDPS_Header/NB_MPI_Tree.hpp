#include <mpi.h>
#include <omp.h>

#include <type_traits>

#include <pvfmm.hpp>

#include "NB_MPI_Node.hpp"

namespace pvfmm {

// the tree type for neighbor operations. must use NB_MPI_Node as node type
// class NB_MPI_Tree : public MPI_Tree<MPI_Node<double>> {

template<class NB_MPI_Node_t>
class NB_MPI_Tree : public MPI_Tree<NB_MPI_Node_t> {
 friend FMM_Pts<typename FMM_Mat_t::FMMNode_t>;

 public:

  typedef typename FMM_Mat_t::FMMNode_t Node_t;
  typedef typename FMM_Mat_t::Real_t Real_t;

  /**
   * \brief Constructor.
   */
  FMM_Tree(MPI_Comm c): MPI_Tree<Node_t>(c), fmm_mat(NULL), bndry(FreeSpace) { };

  /**
   * \brief Virtual destructor.
   */
  virtual ~FMM_Tree(){
  }

  /**
   * \brief Initialize the distributed MPI tree.
   */
  virtual void Initialize(typename Node_t::NodeData* data_) ;

  /**
   * \brief Initialize FMM_Tree.
   */
  void InitFMM_Tree(bool refine, BoundaryType bndry=FreeSpace);

  /**
   * \brief Run FMM
   */
  void SetupFMM(FMM_Mat_t* fmm_mat_);

  /**
   * \brief Run FMM
   */
  void RunFMM();

  /**
   * \brief Clear FMM data: multipole, local expansions and target potential.
   */
  void ClearFMMData();

  /**
   * \brief Build interaction lists for all nodes.
   */
  void BuildInteracLists();

  /**
   * \brief Upward FMM pass (Including MultipoleReduceBcast).
   */
  void UpwardPass();

  /**
   * \brief Reduction and broadcast of multipole expansions.
   */
  void MultipoleReduceBcast() ;

  /**
   * \brief Downward FMM pass.
   */
  void DownwardPass();

  /**
   * \brief Copy FMM output to the tree.
   */
  void Copy_FMMOutput();

 protected:

  std::vector<Matrix<Real_t> > node_data_buff;
  pvfmm::Matrix<Node_t*> node_interac_lst;
  InteracList<Node_t> interac_list;
  FMM_Mat_t* fmm_mat; //Computes all FMM translations.
  BoundaryType bndry;

  std::vector<Matrix<char> > precomp_lst; //Precomputed data for each interaction type.
  std::vector<SetupData<Real_t> > setup_data;

  std::vector<Vector<Real_t> > upwd_check_surf;
  std::vector<Vector<Real_t> > upwd_equiv_surf;
  std::vector<Vector<Real_t> > dnwd_check_surf;
  std::vector<Vector<Real_t> > dnwd_equiv_surf;

};

} // namespace pvfmm