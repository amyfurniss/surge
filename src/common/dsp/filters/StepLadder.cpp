#include "globals.h"
#include "StepLadder.h"
#include "QuadFilterUnit.h"
#include "FilterCoefficientMaker.h"
#include "DebugHelpers.h"
#include "SurgeStorage.h"

namespace StepLadder
{

   void makeCoefficients(FilterCoefficientMaker *cm, float freq, float reso, int sub, SurgeStorage *storage)
   {
      auto pitch = storage->note_to_pitch( freq + 69 ) * Tunings::MIDI_0_FREQ;
      cm->C[sl_cutoff] = fmin(pitch, 22100.0) * 2.0 * M_PI;
      cm->C[sl_reso] = reso * 4; // code says 0-10 is value
   }
   
   __m128 process(QuadFilterUnitState * __restrict f, __m128 inm)
   {
      static constexpr int ssew = 4;
      float in[ssew];
      _mm_store_ps(in, inm);

      float C[2][ssew];
      
      for (int i=0; i < 2; ++i)
      {
         _mm_store_ps(C[i], f->C[i]);
      }
      
      float registers[n_registers][ssew];
      
      for (int i=0; i < 5; ++i)
      {
         _mm_store_ps(registers[i], f->R[i]);
      }
      
       float out[ssew];
       double T = 1.0 / 44100.0;
       for (int v = 0; v < ssew; ++v)
       {
         
         double wa = (2.0 / T) * tan(C[sl_cutoff][v] * T / 2.0);
         double G = wa * T / 2.0;
         double g = G / (1.0 + G);
      
         double u = (in[v] - C[sl_reso][v] * registers[STAGE_4_STORAGE][v]) / (1.0 + g);

         double u1 = (u - registers[STAGE_1_STORAGE][v]) / (1.0 + g);
         double stage_1_out = g * u1 + registers[STAGE_1_STORAGE][v];
         double s1 =  g * u1 + registers[STAGE_1_STORAGE][v];

         double u2 = (stage_1_out - registers[STAGE_2_STORAGE][v]) / (1.0 + g);
         double stage_2_out = g * u2 + registers[STAGE_2_STORAGE][v];
         double s2 =  g * u2 + registers[STAGE_2_STORAGE][v];

         double u3 = (stage_2_out - registers[STAGE_3_STORAGE][v]) / (1.0 + g);
         double stage_3_out = g * u3 + registers[STAGE_3_STORAGE][v];
         double s3 =  g * u3 + registers[STAGE_3_STORAGE][v];

         double u4 = (stage_3_out - registers[STAGE_4_STORAGE][v]) / (1.0 + g);
         double stage_4_out = g * u4 + registers[STAGE_4_STORAGE][v];
         double s4 =  g * u4 + registers[STAGE_4_STORAGE][v];
         out[v] = stage_4_out;

         registers[STAGE_4_OUT][v] = stage_4_out;
         registers[STAGE_1_STORAGE][v] = s1;
         registers[STAGE_2_STORAGE][v] = s2;
         registers[STAGE_3_STORAGE][v] = s3;
         registers[STAGE_4_STORAGE][v] = s4;
      }
   
      __m128 outm = _mm_load_ps(out);
      
      for (int i=0; i<n_registers; ++i)
      {
           f->R[i] = _mm_load_ps(registers[i]);
      }
      return outm;
   }
}
