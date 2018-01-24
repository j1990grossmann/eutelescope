#ifndef DATAFORMAT_TREE_H
#define DATAFORMAT_TREE_H

// c++
#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <vector>

// Root
#include "TObject.h"
#include "TMath.h"

static constexpr unsigned int MEMORY=96;
static constexpr unsigned int CHANNELS=48;
static constexpr unsigned int COLUMNS=16;
static constexpr unsigned int ROWS=3;
static constexpr unsigned int ASSEMBLY=2;
static constexpr unsigned int TIMESTAMP_RANGE=65536;
static constexpr unsigned int THRESHOLD_RANGE=256;

static constexpr double PITCH_STD__X=100;
static constexpr double PITCH_EDGE_X=200;
static constexpr double PITCH_STD_EDGE_Y=1746;
static constexpr double PITCH_STD_Y=1446;
// Weights for COG calculation
static constexpr double WEIGHT_DEFAULT=144600;
static constexpr double WEIGHT_STD_EDGE=289200;
static constexpr double WEIGHT_LOWER_ROW=174600;
static constexpr double WEIGHT_LOWER_CORNER=349200;


// static constexpr float LOWER_ROW=873;
// static constexpr float ROW_1=1596;
// static constexpr float ROW_2=3042;
static constexpr double LOWER_ROW=0.5*PITCH_STD_EDGE_Y;
static constexpr double ROW_1=PITCH_STD_EDGE_Y+0.5*PITCH_STD_Y;
static constexpr double ROW_2=PITCH_STD_EDGE_Y+1.5*PITCH_STD_Y;

static constexpr double TRANSLATE_X=2.*PITCH_EDGE_X+14*PITCH_STD__X;
static constexpr double TRANSLATE_Y=2.*(PITCH_STD_EDGE_Y+2.*PITCH_STD_Y);


// 
struct MemoryNoProcessingBranch_t {
   ULong64_t       pixelMatrix[MEMORY];
   UShort_t        bunchCrossingId[MEMORY];
   UChar_t         header[MEMORY];
   UChar_t         numEvents;
   UChar_t         corrupt;
};


struct RippleCounterBranch_t{
 UInt_t   header;
 UShort_t pixel[CHANNELS];
};
//! my_utilities class
class MemoryCluster: public TObject {
 public:
     double cog_x;
     double cog_y;
     double area;
     unsigned int size_x;
     unsigned int size_y;
     unsigned int clustersize;
     int Chip_Position_Mask;
     int IntraChipPosition;
     unsigned short BX_ID;
     std::set<unsigned int> x_pixel;
     std::set<unsigned int> y_pixel;
  //! my_utilities constructor
  MemoryCluster():
  cog_x(0),
  cog_y(0),
  area(0),
  size_x(0),
  size_y(0),
  clustersize(0),
  Chip_Position_Mask(0),
  IntraChipPosition(0),
  BX_ID(0)
  {};
  ~MemoryCluster(){};
  ClassDef(MemoryCluster,1)
};

class CounterCluster: public TObject {
 public:
     double cog_x;
     double cog_y;
     double area;
     unsigned int size_x;
     unsigned int size_y;
     unsigned int clustersize;
     unsigned int  Chip_Position_Mask;
     int IntraChipPosition;
     std::set<unsigned int> x_pixel;
     std::set<unsigned int> y_pixel;
  //! my_utilities constructor
  CounterCluster():
  cog_x(0),
  cog_y(0),
  area(0),
  size_x(0),
  size_y(0),
  clustersize(0),
  Chip_Position_Mask(0),
  IntraChipPosition(0)
  {};
  ~CounterCluster(){};
  ClassDef(CounterCluster,1)
};
#endif // MY_UTILITIES_H
