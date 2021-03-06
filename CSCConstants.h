#ifndef CSCCONSTANTS_H_
#define CSCCONSTANTS_H_

#include <cmath>
namespace emu { namespace pc {

namespace CSCConstants
{
     enum dcfebs {NUM_DCFEBS = 7};
     enum layer_info  {NUM_LAYERS = 6, KeyLayer = 3};
     enum clctpatterns_info {NUM_PATTERNS = 11};
     enum strips_info {NUM_HSPERDCFEB = 32, NUM_HS_ME1B = 128, NUM_HS_ME1A = 96, NUM_HS = 224, NUM_STRIPS = 112, NUM_DISTRIPS = 56};

int CCLUT_offset[5][6] = {
			{ 0,  0,  0,  0,  0,  0},
			{-2, -1,  0,  0,  1,  2},
			{ 2,  1,  0,  0, -1, -2},
			{-4, -3, -1,  1,  3,  4},
			{ 4,  3,  1, -1, -3, -4}
	};


}
}}
#endif
