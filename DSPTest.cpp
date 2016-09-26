#include <stdio.h>
#include "DSP.h"
#include <math.h>


using namespace DSP;

int main(void){

		samples_t<21,double>  coeficiente;
		samples_t<300> senial, filtrada;
		samples_t<600> asdf;
		
		
		for (int i = -10; i < 11; i++) {
			coeficiente << (10-abs(i))/12.5;
		}
				
		for (int i = 0; i < 300; i++) {
			
			senial << 0.025 * i + 1.25*sin( i*0.75); 
			
		}
		
		filtrada = senial.conv(coeficiente);
		
		return 0;
		
	}

