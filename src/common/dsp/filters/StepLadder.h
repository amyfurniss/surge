#pragma once

class QuadFilterUnitState;
class FilterCoefficientMaker;
class SurgeStorage;

 

namespace StepLadder
{
static constexpr int sl_cutoff = 0, sl_reso = 1;
static int n_registers = 5;


   enum SLRegister {
      STAGE_1_STORAGE = 0,
      STAGE_2_STORAGE = 1,
      STAGE_3_STORAGE = 2,
      STAGE_4_STORAGE = 3,
      STAGE_4_OUT = 4
   };

   void makeCoefficients( FilterCoefficientMaker *cm, float freq, float reso, int subtype, SurgeStorage *storage );
   __m128 process( QuadFilterUnitState * __restrict f, __m128 in );
}
