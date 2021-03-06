// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTelTripletGBLUtility_h
#define EUTelTripletGBLUtility_h 1

#include <memory>
#include "marlin/Processor.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include "lcio.h"
#include <UTIL/CellIDEncoder.h>
#include <UTIL/CellIDDecoder.h>
#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackerHitImpl.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackerDataImpl.h>
#include <IMPL/LCFlagImpl.h>
#include <Exceptions.h>
#include <IMPL/TrackImpl.h>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// ROOT includes
#include <TMatrixD.h>
#include "TH1D.h"

// system includes <>
#include <string>
#include <vector>
#include <map>
#include <deque>

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

// Eigen include
#include <Eigen/Core>
using namespace marlin;

namespace eutelescope {

  // does not need to inherit from marlin::Processor, does it?
  //class EUTelTripletGBLUtility : public marlin::Processor {
  class EUTelTripletGBLUtility {

    public:

      //! Default constructor
      EUTelTripletGBLUtility();

      // Calculate Point-To-Point Jacobian Transport Matrix for distance "ds"
      Eigen::Matrix<double, 5,5> JacobianPointToPoint( double ds );

      // Set your parent
      void setParent(marlin::Processor * par){
	parent = par;
      };

      void bookHistos();

      class hit {
	public:
	  // Coordinates and their position uncertainty
	  double x;
	  double ex;
	  double y;
	  double ey;
	  double z;
	  double ez;
	  // Plane to which the hit belongs
	  unsigned int plane;
	  // clustersize of the cluster associated to the hit
	  int clustersize;
	  // local coords
	  double locx;
	  double locy;
	  int id;
	  // Overloading ostream operator for printing hits:
	  friend std::ostream& operator << (std::ostream& out, const hit& point) // output
	  {
	    out << "(" << point.plane << ", " 
	      << point.x << " +/- " << point.ex << " | " 
	      << point.y << " +/- " << point.ey << " | " 
	      << point.z << " +/- " << point.ez << ")";
	    return out;
	  }
      };

      class triplet {
	public:
	  triplet();
	  triplet(hit hit0, hit hit1, hit hit2);

	  // Keep track of linking status to DUT and REF:
	  bool linked_dut;
	  //bool linked_ref;

	  hit getpoint_at(double z);


	  // Returns x coordinate of the triplet at given z:
	  double getx_at(double z);

	  // Return dx = (x_end - x_start) for the full triplet:
	  double getdx();

	  // Returns dx = (x_measure - x_triplet) in the given plane ipl:
	  double getdx(int ipl);

	  // Returns dx for a given point:
	  double getdx(hit point);


	  // Returns y coordinate of the triplet at given z:
	  double gety_at(double z);

	  // Return dy = (y_end - y_start) for the full triplet:
	  double getdy();

	  // Returns dy = (y_measure - y_triplet) in the given plane ipl:
	  double getdy(int ipl);

	  // Returns dy for a given point:
	  double getdy(hit point);


	  // Return dz = (z_end - z_start) for the full triplet:
	  double getdz();

	  // Returning the hit for the given plane ID
	  hit gethit(int plane);

	  //! Returning the center point of the triplet:
	  hit base();

	  //! Returning the slope of the triplet (x,y):
	  hit slope();

	  friend std::ostream& operator << (std::ostream& out, triplet trip)
	  {
	    out << "Triplet: " << std::endl;
	    for( std::map<unsigned int,hit>::iterator itr = trip.hits.begin(); itr != trip.hits.end(); itr++) {
	      out << "    " << itr->second << std::endl;
	    }
	    return out;
	  };

	private:
	  void filltriplet(hit hit0, hit hit1, hit hit2) {
	    hits.insert( std::pair<unsigned int,hit>(hit0.plane,hit0));
	    hits.insert( std::pair<unsigned int,hit>(hit1.plane,hit1));
	    hits.insert( std::pair<unsigned int,hit>(hit2.plane,hit2));
	  };
	  //! The hits belonging to the triplet:
	  /* Use map since it's already ordered according to plane IDs.
	   * We rely on begin() and rbegin() to deliver pointers to the first and last plane of the triplet.
	   */
	  std::map<unsigned int,hit> hits;   
      };

      class track {
	public:
	  //! Default Track constructor. To be called with two triplets.
	  track(triplet up, triplet down);

	  //! Return the track kink angle in x
	  double kink_x();

	  //! Return the track kink angle in y
	  double kink_y();

	  //! Return the intersection point of the triplets
	  hit intersect();

	  //! Return the track upstream triplet
	  triplet get_upstream();

	  //! Return the track downstream triplet
	  triplet get_downstream();

	  //! Return the track hit in a given plane
	  hit gethit(int plane);

	private:
	  //! Members to store the up- and downstream triplets
	  triplet upstream;
	  triplet downstream;
      };

      //! Find hit triplets from three telescope planes
      /*! This runs over all hits in the planes of the telescope and
       * tries to match triplets by comparing with the middle planes.
       * Two cut criteria can be set:
       * @param triplet_res_cut Cut on the hit residual in the middle plane with respect to the triplet defined by first and last plane
       * @param triplet_slope_cut Cut on the triplet track angle
       *
       * @return a vector of found triplets among the given set of hits.
       */
      void FindTriplets(std::vector<EUTelTripletGBLUtility::hit> &hits, unsigned int plane0, unsigned int plane1, unsigned int plane2, double trip_res_cut, double trip_slope_cut, std::vector<EUTelTripletGBLUtility::triplet> &trip);

      //! Match the upstream and downstream triplets to tracks
      void MatchTriplets(std::vector<EUTelTripletGBLUtility::triplet> &up, std::vector<EUTelTripletGBLUtility::triplet> &down, double z_match, double trip_matching_cut, std::vector<EUTelTripletGBLUtility::track> &track);

      //! Check isolation of triplet within vector of triplets
      bool IsTripletIsolated(std::vector<EUTelTripletGBLUtility::triplet>::iterator it, std::vector<EUTelTripletGBLUtility::triplet> &trip, double z_match, double isolation = 0.3);

      //! Calculate efficiency of plane
      /*! This creates non-standard triplets and driplets (use only 5 planes to contruct them) and looks for matching hit on plane under test
       * Inputs:
       * - triplet
       * - driplet
       * - plane under test (PUT)
       * - z match position
       * - isolation cut
       * - track match cut
       * - Profile that should be filled
       *
       * returns double average efficiency
       */
      double PlaneEfficiency(std::vector<EUTelTripletGBLUtility::triplet> &up, std::vector<EUTelTripletGBLUtility::triplet> &down, std::vector<EUTelTripletGBLUtility::hit> &hits, unsigned int PUT, double track_match_z, double DUT_z, double match_cut, double eff_radius, std::vector<AIDA::IProfile1D*> &profile);

    private:

      //! Fill the telescope plane correlation plots:
      //void TelescopeCorrelationPlots(std::vector<hit> &telescopehits);

      //! Find hit triplets from three telescope planes
      /*! This runs over all hits in the planes of the telescope and
       * tries to match triplets by comparing with the middle planes.
       * Two cut criteria can be set:
       * @param triplet_residual_cut Cut on the hit residual in the middle plane with respect to the triplet defined by first and last plane
       * @param triplet_angle_cut Cut on the triplet track angle
       *
       * @return a vector of found triplets among the given set of hits.
       */
      //void FindTriplets(std::vector<hit> &hits, unsigned int plane0, unsigned int plane1, unsigned int plane2, std::vector<triplet> &trip);

      //! Match the upstream and downstream triplets to tracks
      //void MatchTriplets(std::vector<triplet> &up, std::vector<triplet> &down, double z_match, std::vector<track> &track);

      //! Check isolation of triplet within vector of triplets
      //bool IsTripletIsolated(std::vector<triplet>::iterator it, std::vector<triplet> &trip, double z_match, double isolation = 0.3);

      //! store the parent, needed for having histograms in the same file as the processor that calls the util class
      marlin::Processor * parent;



    protected:
      //! Silicon planes parameters as described in GEAR
      /*! This structure actually contains the following:
       *  @li A reference to the telescope geoemtry and layout
       *  @li An integer number saying if the telescope is w/ or w/o DUT
       *  @li An integer number saying the number of planes in the
       *  telescope.
       *
       *  This object is provided by GEAR during the init() phase and
       *  stored here for local use.
       */
      gear::SiPlanesParameters * _siPlanesParameters;

      //! Silicon plane layer layout
      /*! This is the real geoemetry description. For each layer
       *  composing the telescope the relevant information are
       *  available.
       *
       *  This object is taken from the _siPlanesParameters during the
       *  init() phase and stored for local use
       */
      gear::SiPlanesLayerLayout * _siPlanesLayerLayout;

      std::string _inputCollectionTelescope;

      // Histos

      AIDA::IHistogram1D * sixkxHisto; //driplet-triplet
      AIDA::IHistogram1D * sixkyHisto;
      AIDA::IHistogram1D * sixdxHisto;
      AIDA::IHistogram1D * sixdyHisto;
      AIDA::IHistogram1D * sixdxcHisto;
      AIDA::IHistogram1D * sixdycHisto;

      AIDA::IHistogram1D * sixkxcHisto;
      AIDA::IHistogram1D * sixkycHisto;
      AIDA::IHistogram1D * sixxHisto;
      AIDA::IHistogram1D * sixyHisto;
      AIDA::IHistogram2D * sixxyHisto;
      AIDA::IHistogram2D * sixxycHisto;

      AIDA::IHistogram1D * kinkx;
      AIDA::IHistogram1D * kinky;
      AIDA::IHistogram1D * kinkxy;
      AIDA::IProfile2D * kinkxvsxy;
      AIDA::IProfile2D * kinkyvsxy;
      AIDA::IProfile2D * kinkxyvsxy;

      AIDA::IHistogram1D * triddaMindutHisto;


  };

}

#endif
