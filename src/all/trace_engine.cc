#include "trace_engine.h"
#include "recon_space.h"
#include "sirt.h"
#include "pml.h"
#include "mlem.h"

TraceEngine::TraceEngine(TraceData &trace_data, DISPCommBase<float> &dcomm, TraceRuntimeConfig &conf)
  : trace_data {trace_data},
    comm {dcomm},
    config {conf}
{
  std::string &recon_alg = config.kReconstructionAlg;
  int thread_count = config.thread_count;

  /* SIRT */
  if(recon_alg=="sirt"){
    main_recon_space= new SIRTReconSpace(
        trace_data.metadata().num_slices(), 
        2*trace_data.metadata().num_cols()*trace_data.metadata().num_cols());
    main_recon_space->Initialize(trace_data.metadata().num_grids());
    float init_val=0.;
    main_recon_space->reduction_objects().ResetAllItems(init_val);

    engine.reset(
        new DISPEngineReduction<AReconSpace, float>(
          &comm,
          main_recon_space,
          thread_count));
  }

  /* MLEM */
  else if(recon_alg=="mlem"){
    trace_data.metadata().InitImage(1.);
    main_recon_space= new MLEMReconSpace(
        trace_data.metadata().num_slices(), 
        2*trace_data.metadata().num_cols()*trace_data.metadata().num_cols());
    main_recon_space->Initialize(trace_data.metadata().num_grids());
    float init_val=0.;
    main_recon_space->reduction_objects().ResetAllItems(init_val);

    engine.reset(
        new DISPEngineReduction<AReconSpace, float>(
          &comm,
          main_recon_space,
          thread_count));
  }

  /* PML */
  else if(recon_alg=="pml"){
    trace_data.metadata().InitImage(1.);
    PMLDataRegion *sinograms = new PMLDataRegion(  
        dynamic_cast<DataRegionBase<float, TraceMetadata>&>(
          trace_data.sinograms()
          ));
    trace_data.sinograms(sinograms);
    main_recon_space= new PMLReconSpace(
        trace_data.metadata().num_slices(), 
        2*trace_data.metadata().num_cols()*trace_data.metadata().num_cols());
    main_recon_space->Initialize(trace_data.metadata().num_grids());
    float init_val=0.;
    main_recon_space->reduction_objects().ResetAllItems(init_val);

    /* Prepare processing engine and main reduction space for other threads */
    engine.reset(
        new DISPEngineReduction<AReconSpace, float>(
          &comm,
          main_recon_space,
          thread_count));
  }

  else if(recon_alg=="apmlr"){
    std::cerr << "Algorithm is not ready: " << recon_alg << std::endl;
    exit(0);
  }
  else{
    std::cerr << "Unknown algorithm: " << recon_alg << std::endl;
    exit(0);
  }
}


void TraceEngine::IterativeReconstruction(TraceData &trace_data, int iteration){
  /**************************/
  /* Perform reconstruction */
  /* Define job size per thread request */
  int64_t req_number = trace_data.metadata().num_cols();
  float init_val = 0.;

  auto &main_recon_replica = main_recon_space->reduction_objects();

  #ifdef TIMERON
  std::chrono::duration<double> recon_tot(0.), inplace_tot(0.), update_tot(0.);
  #endif
  for(int i=0; i<iteration; ++i){
    std::cout << "Iteration: " << i << std::endl;
    #ifdef TIMERON
    auto recon_beg = std::chrono::system_clock::now();
    #endif
    engine->RunParallelReduction(trace_data.sinograms(), req_number);  /// Reconstruction
    #ifdef TIMERON
    recon_tot += (std::chrono::system_clock::now()-recon_beg);
    auto inplace_beg = std::chrono::system_clock::now();
    #endif
    engine->ParInPlaceLocalSynchWrapper();              /// Local combination
    #ifdef TIMERON
    inplace_tot += (std::chrono::system_clock::now()-inplace_beg);

    /// Penalized 
    if(config.kReconstructionAlg == "pml"){
      dynamic_cast<PMLDataRegion&>(trace_data.sinograms()).SetFG(0.);
      dynamic_cast<PMLReconSpace&>(*main_recon_space).
        CalculateFG(trace_data.sinograms(), config.b0);
    }

    /// Update reconstruction object
    auto update_beg = std::chrono::system_clock::now();
    #endif
    main_recon_space->UpdateRecon(trace_data, main_recon_replica);
    #ifdef TIMERON
    update_tot += (std::chrono::system_clock::now()-update_beg);
    #endif

    /// Reset iteration
    engine->ResetReductionSpaces(init_val);
    trace_data.sinograms().ResetMirroredRegionIter();
  }
}


void TraceEngine::IterativeReconstruction(){
  IterativeReconstruction(trace_data, config.iteration_count);
}

