// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-
/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */
#ifndef EUTELMULTILINEFIT_H
#define EUTELMULTILINEFIT_H

// built only if GEAR is available
#ifdef USE_GEAR
// eutelescope includes ".h"

// marlin includes ".h"
#include "marlin/Processor.h"

// marlin util includes
#include "Mille.h"
#include "include/MilleBinary.h"

// gear includes <.h>
#include <gear/SiPlanesParameters.h>
#include <gear/SiPlanesLayerLayout.h>

// lcio includes <.h>
#include <EVENT/LCRunHeader.h>
#include <EVENT/LCEvent.h>
#include <IMPL/TrackerHitImpl.h>


#include "EUTelUtility.h"

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#endif

// system includes <>
#include <string>
#include <vector>
#include <map>

#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
#include "TMinuit.h"
#include <TSystem.h>
#include <TMath.h>
class TMinuit;
#endif

#include <iostream>

namespace eutelescope {


//Specify a Rectangular in a sensor
class SensorRectangular {
protected:
    // SensorID
    int sensor;
    // lowest pixel in X direction
    int A;
    // lowest pixel in Y direction
    int B;
    // highest pixel in X direction
    int C;
    // highest pixel in Y direction
    int D;
public:
    SensorRectangular(int s, int a, int b, int c, int d) : sensor(s), A(a), B(b), C(c), D(d) {};
    SensorRectangular() : sensor(0), A(0), B(0), C(0), D(0) {};
    int getSensor() const {
        return sensor;
    }
    //look if x and y are inside the foreseen rectangular
    bool isInside(int x, int y) const {
        return (x >=A && x <=C && y >=B && y <=D);
    }
    void print() {
        streamlog_out(MESSAGE4) << "Sensor: " << sensor << ": (" << A << "|" << B << ") to (" << C << "|" << D << ")" << std::endl;
    }

};

class RectangularArray {
protected:
    std::map<int,SensorRectangular > _rect;

public:
    void addRectangular(SensorRectangular &s) {
        _rect[s.getSensor() ] = s;
    }

    bool isInside(int s, int x, int y) {
        std::map<int,SensorRectangular >::iterator it = _rect.find(s);
        if (it == _rect.end()) { // not in the map means no limit on this sensor -> always true
            return true;
        }
        SensorRectangular cSensor = _rect[s];
        return cSensor.isInside(x,y);
    }

};


class EUTelMilleGBL : public marlin::Processor {

public:
#if defined(USE_ROOT) || defined(MARLIN_USE_ROOT)
    class hit
    {
    public:
        hit()
        {
        }
        hit(double tx, double ty, double tz, double rx, double ry, double rz,int i)
        {
            x = tx;
            y = ty;
            z = tz;
            resolution_x = rx;
            resolution_y = ry;
            resolution_z = rz;

            planenumber = i;
        }
        double x;
        double y;
        double z;
        double resolution_x;
        double resolution_y;
        double resolution_z;

        int planenumber;
    };
#endif

    //! Variables for hit parameters
    class HitsInPlane {
    public:
        HitsInPlane() {
            measuredX = 0.0;
            measuredY = 0.0;
            measuredZ = 0.0;
        }
        HitsInPlane(double x, double y, double z)
        {
            measuredX = x;
            measuredY = y;
            measuredZ = z;
        }
        bool operator<(const HitsInPlane& b) const
        {
            return (measuredZ < b.measuredZ);
        }
        double measuredX;
        double measuredY;
        double measuredZ;
    };
//     Template trajectory helper classes
//     Measurement or Scatterer at GBL Point
    enum class GBL_Type : uint8_t {
        scatterer,
        measurement,
        air
    };
    friend std::ostream& operator<< (std::ostream& rStream, const GBL_Type& flag ) {
        switch (flag) {
        case GBL_Type::scatterer :
            rStream << "scatterer";
            break;
        case GBL_Type::measurement :
            rStream << "measurement";
            break;
        case GBL_Type::air :
            rStream << "air";
            break;
        default :
            rStream << "Not a GBL_TYPE";
            break;
        }
        return rStream;
    }
//     Measurement included in the Fit Biased/Unbiased residuals (See GBL description)
    enum class GBL_Active : uint8_t {
        active,
        not_active
    };
    friend std::ostream& operator<< (std::ostream& rStream, const GBL_Active& flag ) {
        switch (flag) {
        case GBL_Active::active :
            rStream << "active";
            break;
        case GBL_Active::not_active :
            rStream << "not active";
            break;
        default :
            rStream << "Not a GBL_Active flag";
            break;
        }
        return rStream;
    }
    class GBL_Point_Template {
        //         GBL_Point_Template(){
        //         }
    public:
        double _step;
        double _InverseScatAngleSqrt;
        double _zpos;
        double _radlen;
        GBL_Type   _gbl_type;
        GBL_Active _gbl_active;
    };
    
//  Object holding information about scattering Material to construct a GBL trajectory with a material budget estimate
//  Interfaces with the gear file geometry and is used to calculate the Scattering material
    class GBL_Trajectory_Template {
    public:
        GBL_Trajectory_Template():
            _distplane(0), _kappa(1.), _p(0), _sumeps(0), _n_dut(0), _n_deadlayer(0), _n_telplanes(0)
        {
            _proL2m = Eigen::Matrix2d::Identity();
            //             Zero initialize measurement precision 1/resolution^2
            _measPrec = Eigen::Vector2d::Zero();
            //             Mean scattering angle is zero
            _scat = Eigen::Vector2d::Zero();

        }
        //! Prints the tracker geometry from a gear file.
        /*! This prints the tracker geometry.
         */
        void print_geometry(const gear::SiPlanesLayerLayout&  siplaneslayerlayout);
        //! Prints the template trajectory material budged and inverse scattering angles as calculated from gear file.
        /*! This prints the track template used to fill a GBL trajectory
         */
        void print_template_traj();
        //! Called for every run.
        /*! This is executed at the beginning of every run and fills a object containing the static gemoetry
         */
        void fill_geometry(const gear::SiPlanesLayerLayout&  siplaneslayerlayout);
        //! Called for every run.
        /*! Fill measurement precision vector from previous iteration with values obtained with Millepede
         */
        void set_res(double res_x, double res_y);
        //! Called for every run.
        /*! Fill measurement precision vector. Telescope resolution estimation with momentum in GeV and alignement steps.
         */
        void set_res(double p_momentum, int first_alignment_step);
        //! Called for every run.
        /*! Set kappa for highland formula default is 1
         */
        void set_kappa(double kappa);
        void set_p(double p);
        double get_kappa();
        double get_p();
        const std::vector<GBL_Point_Template>& getTemplateTraj() const { return _GBL_Template_Vec; }
        const Eigen::Matrix2d& getproL2m() const { return _proL2m; }
        const Eigen::Vector2d& getmeasPrec() const { return _measPrec; }
        const Eigen::Vector2d& getscat() const { return _scat; }
        //! Called for every run.
        /*! This is executed for every track fed into the GBL fit routine.
         *  Fills an object with the absorber material and  central part of the
         *  angle distribution according to the gemoetry given in a gear file and a given track candidate with the highland formula.
         *  FIXME The material thickness accounts for the track angle with respect to the plane of absorber material.
         */
        //! Called for every run.
        /*! This is executed at the beginning of every run and fills a object containing the static gemoetry
         */
        //         void fill_materialbudget(gear::SiPlanesLayerLayout * const siplaneslayerlayout);
        /*! Prints the material budged and the GBL trajectory for a given track candidate.
         */
//     protected:
    private:
        void fill_sort_ids(const gear::SiPlanesLayerLayout& siplaneslayerlayout);
        void fill_x0(const gear::SiPlanesLayerLayout& siplaneslayerlayout, const std::vector<double>& layer_dist);
//         See paper "Performance of the EUDET-type beam telescopes" equation 2
        inline double scattering_fraction_var(double eps_scat);
//      Refer to the first measurement or scatter z position upstream and the last Measurement downstream
        double _distplane;
        double _kappa;
        double _p;
        double _sumeps ;
        unsigned int _n_dut, _n_deadlayer, _n_telplanes;
        std::vector< std::pair<unsigned int, double>> _id_zpos_vec;
//         This vector stores the data used to build up the gbl trajectory
        std::vector< GBL_Point_Template > _GBL_Template_Vec;

        Eigen::Matrix2d _proL2m;
        Eigen::Vector2d _measPrec;
        Eigen::Vector2d _scat;
        static constexpr double _epsAir = 303900.;
        static constexpr double _airscat_dist_hi  = .5+std::pow(12,-.5);
        static constexpr double _airscat_dist_lo  = .5-std::pow(12,-.5);
    };


    //! Returns a new instance of EUTelMilleGBL
    /*! This method returns a new instance of this processor.  It is
     *  called by Marlin execution framework and it shouldn't be
     *  called/used by the final user.
     *
     *  @return a new EUTelMilleGBL.
     */
    virtual Processor * newProcessor() {
        return new EUTelMilleGBL;
    }

    //DP virtual bool hitContainsHotPixels( TrackerHitImpl   * hit) ;


    //! Default constructor
    EUTelMilleGBL ();

    //! Called at the job beginning.
    /*! This is executed only once in the whole execution. It prints
     *  out the processor parameters and check that the GEAR
     *  environment is properly set up and accessible from Marlin.
     */
    virtual void init ();

    //! Called for every run.
    /*! It is called for every run, and consequently the run counter
     *  is incremented. The geometry ID of the file is compared with
     *  the one provided by the GEAR geometry description. In case the
     *  two are different, the user is asked to decide to quit or to
     *  continue with a description that might be wrong.
     *
     *  @param run the LCRunHeader of the this current run
     */
    virtual void processRunHeader (LCRunHeader * run);

    //! Called for first event per run
    /*! Reads hotpixel information from hotPixelCollection into hotPixelMap
     * to be used in the sensor exclusion area logic
     */
    //DP virtual void  FillHotPixelMap(LCEvent *event);

    //! Called every event
    /*! This is called for each event in the file. Each element of the
     *  pulse collection is scanned and the center of the cluster is
     *  translated into the external frame of reference thanks to the
     *  GEAR geometry description.
     *
     *  The cluster center might be calculate using a standard linear
     *  charge center of gravity algortihm or applying a more
     *  sophisticated non linear eta function. This behaviour is
     *  regulated by the user from the steering file.
     *
     *  @throw UnknownDataTypeException if the cluster type is unknown
     *
     *  @param evt the current LCEvent event as passed by the
     *  ProcessMgr
     */
    virtual void processEvent (LCEvent * evt);


    //! Called after data processing.
    /*! This method is called when the loop on events is
     *  finished.
     */
    virtual void end();


    //! Histogram booking
    /*! Some control histograms are filled during this procedure in
     *  order to be able to perform easy check on the quality of the
     *  output hits and also to understand if the frame of reference
     *  conversion has been properly done. Of course this method is
     *  effectively doing something only in the case MARLIN_USE_AIDA.
     */
    void bookHistos();


protected:


    //! Ordered sensor ID
    /*! Within the processor all the loops are done up to _nPlanes and
     *  according to their position along the Z axis (beam axis).
     *
     *  This vector is containing the sensorID sorted according to the
     *  same rule.
     */
    std::vector< int > _orderedSensorID;
    std::vector< int > _orderedSensorID_wo_excluded;



    //! TrackerHit collection name
    /*! Input collection with hits.
     */
    std::vector<std::string > _hitCollectionName;

    //! TRACK collection name
    /*! Output collection with fitted tracks.
     */
    std::string _trackCollectionName;

    //! Hot pixel collection name.
    /*!
     * this collection is saved in a db file to be used at the clustering level
     */
    //DP std::string _hotPixelCollectionName;

    //! Vector of map arrays, keeps record of hit pixels
    /*! The vector elements are sorted by Detector ID
     *  For each Detector unique ID element a map of pixels is created.
     *  first level key   sensor unique
     *              value sensor map
     *  sensor map key    unique row number
     *             value  vector of column numbers.
     */

    //DP std::map<std::string, bool > _hotPixelMap;


    // parameters

    float _distanceMax;
    std::vector<float> _distanceMaxVec;
    std::vector<int > _excludePlanes; //only for internal usage
    std::vector<int > _excludePlanes_sensorIDs; //this is going to be
    //set by the user.
    std::vector<int > _FixedPlanes; //only for internal usage
    std::vector<int > _FixedPlanes_sensorIDs; //this is going to be
    //set by the user.

    double _eBeam; // DP
    double _triCut;
    double _driCut;
    double _sixCut;
    double _slopeCut;
    double _chi2Cut;
    int _sensorID_DUT;

    int _IsFirstAlignStep;

    double _kappa;
    double _targetthick;

    int _maxTrackCandidates;
    int _maxTrackCandidatesTotal;

    std::string _binaryFilename;

    float _telescopeResolution;
    int _onlySingleHitEvents;
    int _onlySingleTrackEvents;
    Utility::alignMode _alignMode;
    std::string _alignModeString;
    int _useResidualCuts;

    std::vector<float > _residualsXMin;
    std::vector<float > _residualsYMin;
    std::vector<float > _residualsXMax;
    std::vector<float > _residualsYMax;

    std::vector<float> _resolutionX;
    std::vector<float> _resolutionY;
    std::vector<float> _resolutionZ;

    std::vector<int> _sensorID_excluded;
    std::vector<int> _FixParameter;


    int _generatePedeSteerfile;
    std::string _pedeSteerfileName;
    int _runPede;
    int _usePedeUserStartValues;
    std::vector<float > _pedeUserStartValuesX;
    std::vector<float > _pedeUserStartValuesY;
    std::vector<float > _pedeUserStartValuesZ;

    std::vector<float > _pedeUserStartValuesAlpha;
    std::vector<float > _pedeUserStartValuesBeta;
    std::vector<float > _pedeUserStartValuesGamma;

    int _inputMode;
    float _testModeSensorResolution;
    float _testModeXTrackSlope;
    float _testModeYTrackSlope;

    std::vector<float > _testModeSensorZPositions;

    std::vector<float > _testModeSensorXShifts;
    std::vector<float > _testModeSensorYShifts;
    std::vector<float > _testModeSensorGamma;
    std::vector<float > _testModeSensorAlpha;
    std::vector<float > _testModeSensorBeta;

    std::string _alignmentConstantLCIOFile;
    std::string _alignmentConstantCollectionName;

    std::vector<int> _useSensorRectangular;

private:

    //! Run number
    int _iRun;

    //! Event number
    int _iEvt;
    unsigned int _nTri;
    unsigned int _nDri;
    unsigned int _nSix;

    //! counter for printed events (for debugging)
    int _printEventCounter;

    // n Tscope planes
    int _nTelPlanes;

    // Excluded planes
    int _nExcludePlanes;

    // Statistics
    int _nMilleDataPoints;
    int _nMilleTracks;

    // Mille
    //DP Mille * _mille;

    //! Conversion ID map.
    /*! In the data file, each cluster is tagged with a detector ID
     *  identify the sensor it belongs to. In the geometry
     *  description, there are along with the sensors also "passive"
     *  layers and other stuff. Those are identify by a layerindex. So
     *  we need a conversion table to go from the detectorID to the
     *  layerindex.
     */
    std::map< int, int > _conversionIdMap;

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

    // Partly outdated GEAR readings:
    int * _planeSort;
    int * _planeID;
    double * _planePosition;
    double * _planeThickness;
    double * _planeX0;
    double * _planeResolution;

#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
    //! AIDA histogram map
    /*! Instead of putting several pointers to AIDA histograms as
     *  class members, histograms are booked in the init() method and
     *  their pointers are inserted into this map keyed by their
     *  names.
     *  The histogram filling can proceed recalling an object through
     *  its name
     */
    std::map<std::string, AIDA::IBaseHistogram * > _aidaHistoMap;

    static std::string _numberTracksLocalname;

    static std::string _chi2XLocalname;
    static std::string _chi2YLocalname;

    static std::string _residualXLocalname;
    static std::string _residualYLocalname;
    static std::string _residualZLocalname;
#endif

    int _nPlanes;

    std::vector<std::vector<double> > _xPos;
    std::vector<std::vector<double> > _yPos;
    std::vector<std::vector<double> > _zPos;

    double * _xPosHere;
    double * _yPosHere;
    double * _zPosHere;
    double * _waferResidX;
    double * _waferResidY;
    double * _waferResidZ;
    double * _telescopeResolX;
    double * _telescopeResolY;
    double * _telescopeResolZ;
    double * _xFitPos;
    double * _yFitPos;

    std::vector<double> _siPlaneZPosition;
    
    GBL_Trajectory_Template _gbltemplate;
    
    //! Fill histogram switch
    /*! Only for debug reason
     */
    bool _histogramSwitch;

    //! Limits the pixels on each sensor-plane to a sub-rectangular
    RectangularArray _rect;

};

//! A global instance of the processor
EUTelMilleGBL gEUTelMilleGBL;

}
#endif
#endif
