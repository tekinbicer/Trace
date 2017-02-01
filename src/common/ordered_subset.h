#include <vector>
#include <cmath>
#include <cstring>
#include <stdexcept>
#include "trace_data.h"
#include "data_region_base.h"
#include "disp_engine_reduction.h"

enum class SubsetType {Single, Ordered, Distributed};

class DataSubset
{
  private:
    float *projs_;
    float *thetas_;
    int nprojs_; 
    int ntotsinogs_;
    int ncols_;
    int nsubsets_;
    int nsinogs_;
    int beg_sinog_id_;
    int center_;
    int nneighbors_;
    int recon_init_val_;
    size_t ray_per_proj_;

    std::vector<float*> subsets_;
    std::vector<float*> subset_thetas_;
    std::vector<std::vector<int>> indices_;

    DataRegionBareBase<float> *recon_;

    void GenerateOrderedIndices() {
      std::vector<int> p;

      for(int i=0; i<nprojs_; ++i) 
        p.push_back(i);

      int m=nsubsets_;
      for(int i=0; i<m; ++i){
        std::vector<int> k {};
        indices_.push_back(k);
      }

      for(int i=0; i<p.size(); ++i)
        indices_[i%m].push_back(p[i]);
    }

    void OrganizeOrderedSubsets()
    {
      for(int i=0; i<nsubsets_; ++i){
        std::cout << "allocating=" << ray_per_proj_ << "x" << indices_[i].size() << 
          "=" << ray_per_proj_*indices_[i].size() << " fp # for subset=" << i << '\n'; 
        subsets_.push_back(new float[ray_per_proj_*indices_[i].size()]);
        subset_thetas_.push_back(new float[indices_[i].size()]);
      }

      for(int i=0; i<indices_.size(); ++i){
        for(int j=0; j<indices_[i].size(); ++j){
          float *source = projs_+indices_[i][j]*ray_per_proj_;
          float *dest = subsets_[i]+j*ray_per_proj_;
          std::cout << "Copying proj=" << indices_[i][j] << " to dest=" << subsets_[i] << "+" << j << "*" << ray_per_proj_  << "=" << dest << std::endl;

          for(size_t k=0; k<ray_per_proj_; ++k)
            dest[k] = source[k];
          
          subset_thetas_[i][j] = thetas_[indices_[i][j]];

          std::cout << "Indices[" << i << "]["<<j<<"]="<<indices_[i][j];
          std::cout << "; Theta[" << i << "]["<<j<<"]="<<subset_thetas_[i][j] << std::endl;
        }
      }
    }

    void ReOrganizeData(SubsetType type) {
      switch(type){
        case SubsetType::Single:
          std::cout << "single subset" << std::endl;
          GenerateOrderedIndices();
          OrganizeOrderedSubsets();
          break;
        case SubsetType::Ordered: 
          // Regular data organization
          std::cout << "ordered subset" << std::endl;
          GenerateOrderedIndices();
          OrganizeOrderedSubsets();
          break;
        case SubsetType::Distributed:
          std::cout << "not implemented yet" << std::endl;
          // Distributed data organization
          //GenerateDistributedIndices();
          //OrganizeDistributedSubsets();
          break;
        default:
          std::cout << "unknown subset type" << std::endl;
          break;
      }
    }

  public:

    DataSubset(
        float *projections, 
        float *thetas,
        int *dims,
        int nsinogs,
        int beg_sinog_id,
        float center): 
      DataSubset(
          projections, thetas, dims, nsinogs, beg_sinog_id, 
          1, center, 0, 1., SubsetType::Single) {}

    DataSubset(
        float *projections, 
        float *thetas,
        int *dims,
        int nsinogs,
        int beg_sinog_id,
        int nsubsets,
        float center): 
      DataSubset(
          projections, thetas, dims, nsinogs, beg_sinog_id, 
          nsubsets, center, 0, 1., SubsetType::Ordered) {}

    DataSubset(
        float *projections, 
        float *thetas,
        int *dims,
        int nsinogs,
        int beg_sinog_id,
        int nsubsets,
        float center,
        float recon_init_val): 
      DataSubset(
          projections, thetas, dims, nsinogs, beg_sinog_id, 
          nsubsets, center, 0, recon_init_val, SubsetType::Ordered) {}

    DataSubset(
        float *projections, 
        float *thetas,
        int *dims,
        int nsinogs,
        int beg_sinog_id,
        int nsubsets,
        float center,
        int nneighbors,
        float recon_init_val,
        SubsetType type
        ):
      projs_ {projections},
      thetas_ {thetas},
      nprojs_ {dims[0]},
      ntotsinogs_ {dims[1]},
      ncols_ {dims[2]},
      nsinogs_ {nsinogs},
      beg_sinog_id_ {beg_sinog_id},
      nsubsets_ {nsubsets},
      center_ {center},
      nneighbors_ {nneighbors},
      recon_init_val_ {recon_init_val},
      ray_per_proj_ {nsinogs_*ncols_}
    {
      // Setup common reconstruction object
      int num_recon_slices = nsinogs_+2*nneighbors_;
      recon_ = new DataRegionBareBase<float>(num_recon_slices*ncols_*ncols_);
      for(int i=0; i<num_recon_slices*ncols_*ncols_; ++i)
        (*recon_)[i] = recon_init_val_;

      std::cout << "inited recon=" << &((*recon_)[0]) << std::endl;

      ReOrganizeData(type);
    }

    ~DataSubset() {
      delete recon_;
      for(int i=0; i<subsets_.size(); ++i){
        delete[] subsets_[i]; 
        delete[] subset_thetas_[i];
      }
    }

    std::vector<float*> Subsets() const { return subsets_; }
    float* Subsets(int id) const { return subsets_[id]; }
    std::vector<float*> SubsetThetas() const { return subset_thetas_; }
    float* SubsetThetas(int id) const { return subset_thetas_[id]; }
    std::vector<std::vector<int>> Indices() const { return indices_; }
    std::vector<int> Indices(int id) const { return indices_[id]; }
    int NSubsets() const { return nsubsets_; }

    DataRegionBase<float, TraceMetadata>* DataRegionSubset(int id) {
      if(id>=nsubsets_)
        throw std::out_of_range(
            "id>=nsubsets: id=" + std::to_string(id) + 
            "; nsubsets=" + std::to_string(nsubsets_));

      TraceMetadata *metadata = new TraceMetadata(
          SubsetThetas(id),
          0, beg_sinog_id_, 0,
          ntotsinogs_, indices_[id].size(), nsinogs_, ncols_, ncols_, 
          center_, nneighbors_, recon_init_val_);
      metadata->recon(*recon_);  // Forward reconstruction object

      std::cout << "returning subset=" << id << "; offset=" << Subsets(id) << "; count=" << metadata->count() << "; recon=" << &((*recon_)[0]) << std::endl;

      DataRegionBase<float, TraceMetadata> *subset_region = 
        new DataRegionBase<float, TraceMetadata>(
            Subsets(id), 
            metadata->count(), metadata);

      return subset_region;
    }
};
