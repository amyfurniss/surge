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
      cm->C[sl_cutoff] = tan(fmin(pitch, 22100.0) * M_PI * dsamplerate_inv);
      cm->C[sl_reso] = reso * 4;
      cm->C[sl_poles] = sub;
   }

   __m128 process(QuadFilterUnitState * __restrict f, __m128 inm)
   {
     
   float C[4];
   _mm_store_ps(C, f->C[sl_poles]);
   // assume that all voices have the same number of poles
   switch ((int)(C[0]))
   {
      case 0:
         return process_1_pole(f, inm);
      case 1:
         return process_2_poles(f, inm);
      case 2:
         return process_3_poles(f, inm);
      case 3:
      default:
         return process_4_poles(f, inm);
     }
   }

   __m128 process_4_poles(QuadFilterUnitState * __restrict f, __m128 inm)
   {
      const __m128 T = _mm_load_ps(new float[]{(float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv});
      const __m128 one = _mm_load_ps(new float[]{1.0f, 1.0f, 1.0f, 1.0f});
      const __m128 two = _mm_load_ps(new float[]{2.0f, 2.0f, 2.0f, 2.0f});
      const __m128 two_div_T =_mm_div_ps(two, T);
      //double wa = (2.0 / T) * tan(C[sl_cutoff] * T / 2.0);
      __m128 wa =_mm_mul_ps(two_div_T, f->C[sl_cutoff]);
      //double G = wa * T / 2.0;
      __m128 G =_mm_div_ps(_mm_mul_ps(wa, T), two);
      //double g = G / (1.0 + G);
      __m128 g =_mm_div_ps(G, _mm_add_ps(one, G));
      //  (1.0 + g)
      __m128 one_plus_g =_mm_add_ps(one, g);

      // u = (inm - C[sl_reso][v] * registers[STAGE_4_STORAGE][v]) / (1.0 + g);
      __m128 u =_mm_div_ps(_mm_sub_ps(inm, _mm_mul_ps(f->C[sl_reso], f->R[LAST_OUTPUT])), one_plus_g);
      //double u1 = (u - registers[STAGE_1_STORAGE][v]) / (1.0 + g);
      __m128 u1 =_mm_div_ps(_mm_sub_ps(u, f->R[STAGE_1_STORAGE]), one_plus_g);
      //double stage_1_out = g * u1 + registers[STAGE_1_STORAGE][v];
      __m128 stage_1_out =_mm_add_ps(_mm_mul_ps(g, u1), f->R[STAGE_1_STORAGE]);
      __m128 s1 =  stage_1_out;

      //double u2 = (stage_1_out - registers[STAGE_2_STORAGE][v]) / (1.0 + g);
      __m128 u2 =_mm_div_ps(_mm_sub_ps(stage_1_out, f->R[STAGE_2_STORAGE]), one_plus_g);
      //double stage_2_out = g * u2 + registers[STAGE_2_STORAGE][v];
      __m128 stage_2_out =_mm_add_ps(_mm_mul_ps(g, u2), f->R[STAGE_2_STORAGE]);
      __m128 s2 =  stage_2_out;
       
      //double u3 = (stage_2_out - registers[STAGE_3_STORAGE][v]) / (1.0 + g);
      __m128 u3 =_mm_div_ps(_mm_sub_ps(stage_2_out, f->R[STAGE_3_STORAGE]), one_plus_g);
      //double stage_3_out = g * u3 + registers[STAGE_3_STORAGE][v];
      __m128 stage_3_out =_mm_add_ps(_mm_mul_ps(g, u3), f->R[STAGE_3_STORAGE]);
      __m128 s3 =  stage_3_out;

      //double u4 = (stage_3_out - registers[STAGE_4_STORAGE][v]) / (1.0 + g);
      __m128 u4 =_mm_div_ps(_mm_sub_ps(stage_3_out, f->R[STAGE_4_STORAGE]), one_plus_g);
      //double stage_4_out = g * u4 + registers[STAGE_4_STORAGE][v];
      __m128 stage_4_out =_mm_add_ps(_mm_mul_ps(g, u4), f->R[STAGE_4_STORAGE]);
      __m128 s4 =  stage_4_out;
      __m128 outm = stage_4_out;

      f->R[LAST_OUTPUT] = stage_4_out;
      f->R[STAGE_1_STORAGE] = s1;
      f->R[STAGE_2_STORAGE] = s2;
      f->R[STAGE_3_STORAGE] = s3;
      f->R[STAGE_4_STORAGE] = s4;
      
      return outm;
   }


   __m128 process_3_poles(QuadFilterUnitState * __restrict f, __m128 inm)
   {
      const __m128 T = _mm_load_ps(new float[]{(float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv});
      const __m128 one = _mm_load_ps(new float[]{1.0f, 1.0f, 1.0f, 1.0f});
      const __m128 two = _mm_load_ps(new float[]{2.0f, 2.0f, 2.0f, 2.0f});
      const __m128 two_div_T =_mm_div_ps(two, T);
      //double wa = (2.0 / T) * tan(C[sl_cutoff] * T / 2.0);
      __m128 wa =_mm_mul_ps(two_div_T, f->C[sl_cutoff]);
      //double G = wa * T / 2.0;
      __m128 G =_mm_div_ps(_mm_mul_ps(wa, T), two);
      //double g = G / (1.0 + G);
      __m128 g =_mm_div_ps(G, _mm_add_ps(one, G));
      //  (1.0 + g)
      __m128 one_plus_g =_mm_add_ps(one, g);

      // u = (inm - C[sl_reso][v] * registers[STAGE_4_STORAGE][v]) / (1.0 + g);
      __m128 u =_mm_div_ps(_mm_sub_ps(inm, _mm_mul_ps(f->C[sl_reso], f->R[LAST_OUTPUT])), one_plus_g);
      //double u1 = (u - registers[STAGE_1_STORAGE][v]) / (1.0 + g);
      __m128 u1 =_mm_div_ps(_mm_sub_ps(u, f->R[STAGE_1_STORAGE]), one_plus_g);
      //double stage_1_out = g * u1 + registers[STAGE_1_STORAGE][v];
      __m128 stage_1_out =_mm_add_ps(_mm_mul_ps(g, u1), f->R[STAGE_1_STORAGE]);
      __m128 s1 =  stage_1_out;

      //double u2 = (stage_1_out - registers[STAGE_2_STORAGE][v]) / (1.0 + g);
      __m128 u2 =_mm_div_ps(_mm_sub_ps(stage_1_out, f->R[STAGE_2_STORAGE]), one_plus_g);
      //double stage_2_out = g * u2 + registers[STAGE_2_STORAGE][v];
      __m128 stage_2_out =_mm_add_ps(_mm_mul_ps(g, u2), f->R[STAGE_2_STORAGE]);
      __m128 s2 =  stage_2_out;
       
      //double u3 = (stage_2_out - registers[STAGE_3_STORAGE][v]) / (1.0 + g);
      __m128 u3 =_mm_div_ps(_mm_sub_ps(stage_2_out, f->R[STAGE_3_STORAGE]), one_plus_g);
      //double stage_3_out = g * u3 + registers[STAGE_3_STORAGE][v];
      __m128 stage_3_out =_mm_add_ps(_mm_mul_ps(g, u3), f->R[STAGE_3_STORAGE]);
      __m128 s3 =  stage_3_out;

      __m128 outm = stage_3_out;

      f->R[LAST_OUTPUT] = stage_3_out;
      f->R[STAGE_1_STORAGE] = s1;
      f->R[STAGE_2_STORAGE] = s2;
      f->R[STAGE_3_STORAGE] = s3;
      
      return outm;
   }


   __m128 process_2_poles(QuadFilterUnitState * __restrict f, __m128 inm)
   {
      const __m128 T = _mm_load_ps(new float[]{(float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv});
      const __m128 one = _mm_load_ps(new float[]{1.0f, 1.0f, 1.0f, 1.0f});
      const __m128 two = _mm_load_ps(new float[]{2.0f, 2.0f, 2.0f, 2.0f});
      const __m128 two_div_T =_mm_div_ps(two, T);
      //double wa = (2.0 / T) * tan(C[sl_cutoff] * T / 2.0);
      __m128 wa =_mm_mul_ps(two_div_T, f->C[sl_cutoff]);
      //double G = wa * T / 2.0;
      __m128 G =_mm_div_ps(_mm_mul_ps(wa, T), two);
      //double g = G / (1.0 + G);
      __m128 g =_mm_div_ps(G, _mm_add_ps(one, G));
      //  (1.0 + g)
      __m128 one_plus_g =_mm_add_ps(one, g);

      // u = (inm - C[sl_reso][v] * registers[STAGE_4_STORAGE][v]) / (1.0 + g);
      __m128 u =_mm_div_ps(_mm_sub_ps(inm, _mm_mul_ps(f->C[sl_reso], f->R[LAST_OUTPUT])), one_plus_g);
      //double u1 = (u - registers[STAGE_1_STORAGE][v]) / (1.0 + g);
      __m128 u1 =_mm_div_ps(_mm_sub_ps(u, f->R[STAGE_1_STORAGE]), one_plus_g);
      //double stage_1_out = g * u1 + registers[STAGE_1_STORAGE][v];
      __m128 stage_1_out =_mm_add_ps(_mm_mul_ps(g, u1), f->R[STAGE_1_STORAGE]);
      __m128 s1 =  stage_1_out;

      //double u2 = (stage_1_out - registers[STAGE_2_STORAGE][v]) / (1.0 + g);
      __m128 u2 =_mm_div_ps(_mm_sub_ps(stage_1_out, f->R[STAGE_2_STORAGE]), one_plus_g);
      //double stage_2_out = g * u2 + registers[STAGE_2_STORAGE][v];
      __m128 stage_2_out =_mm_add_ps(_mm_mul_ps(g, u2), f->R[STAGE_2_STORAGE]);
      __m128 s2 =  stage_2_out;
       
      __m128 outm = stage_2_out;

      f->R[LAST_OUTPUT] = stage_2_out;
      f->R[STAGE_1_STORAGE] = s1;
      f->R[STAGE_2_STORAGE] = s2;
      
      return outm;
   }

   __m128 process_1_pole(QuadFilterUnitState * __restrict f, __m128 inm)
   {
      const __m128 T = _mm_load_ps(new float[]{(float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv, (float)dsamplerate_inv});
      const __m128 one = _mm_load_ps(new float[]{1.0f, 1.0f, 1.0f, 1.0f});
      const __m128 two = _mm_load_ps(new float[]{2.0f, 2.0f, 2.0f, 2.0f});
      const __m128 two_div_T =_mm_div_ps(two, T);
      //double wa = (2.0 / T) * tan(C[sl_cutoff] * T / 2.0);
      __m128 wa =_mm_mul_ps(two_div_T, f->C[sl_cutoff]);
      //double G = wa * T / 2.0;
      __m128 G =_mm_div_ps(_mm_mul_ps(wa, T), two);
      //double g = G / (1.0 + G);
      __m128 g =_mm_div_ps(G, _mm_add_ps(one, G));
      //  (1.0 + g)
      __m128 one_plus_g =_mm_add_ps(one, g);

      //double u1 = (u - registers[STAGE_1_STORAGE][v]) / (1.0 + g);
      __m128 u1 =_mm_div_ps(_mm_sub_ps(inm, f->R[LAST_OUTPUT]), one_plus_g);
      //double stage_1_out = g * u1 + registers[STAGE_1_STORAGE][v];
      __m128 stage_1_out =_mm_add_ps(_mm_mul_ps(g, u1), f->R[STAGE_1_STORAGE]);
      __m128 s1 =  stage_1_out;

      __m128 outm = stage_1_out;

      f->R[LAST_OUTPUT] = stage_1_out;
      f->R[STAGE_1_STORAGE] = s1;
      return outm;
   }

}
