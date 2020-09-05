#pragma once

class QuadFilterUnitState;
class FilterCoefficientMaker;
class SurgeStorage;

 

namespace StepLadder
{
static constexpr int sl_cutoff = 0, sl_reso = 1, sl_poles = 2;
static int n_registers = 5;


   enum SLRegister {
      STAGE_1_STORAGE = 0,
      STAGE_2_STORAGE = 1,
      STAGE_3_STORAGE = 2,
      STAGE_4_STORAGE = 3,
      LAST_OUTPUT = 4
   };

   void makeCoefficients( FilterCoefficientMaker *cm, float freq, float reso, int subtype, SurgeStorage *storage );
   __m128 process( QuadFilterUnitState * __restrict f, __m128 in );
__m128 process_1_pole( QuadFilterUnitState * __restrict f, __m128 in );
__m128 process_2_poles( QuadFilterUnitState * __restrict f, __m128 in );
__m128 process_3_poles( QuadFilterUnitState * __restrict f, __m128 in );
__m128 process_4_poles( QuadFilterUnitState * __restrict f, __m128 in );
}
