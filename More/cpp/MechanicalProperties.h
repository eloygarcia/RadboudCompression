#ifndef __mechanicalProperties_h
#define __mechanicalProperties_h

float bulkModulus( float Eyoung, float nuPoisson){ return (Eyoung/(3*(1-2*nuPoisson))); } ;

float shearModulus( float Eyoung, float nuPoisson){ return (Eyoung/(2*(1+nuPoisson))); } ;

float firstLameParameter( float Eyoung, float nuPoisson){ return ((Eyoung*nuPoisson)/((1+nuPoisson)*(1-2*nuPoisson))); } ;

float secondLameParameter( float Eyoung, float nuPoisson){ return (Eyoung/(2*(1+nuPoisson))); } ;

float anisotropyTerm( float Eyoung, float anisotropyRatio) { return (Eyoung*anisotropyRatio); };


#endif