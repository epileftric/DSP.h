# DSP.h

This namespace defines and implements the template  samples_t<n,type> for signal processing.
The general idea of it's implementation is to simplify the process of dealing with samples
and filtering. 
  
Just by setting two samples_t you can use one for the samples you read and another for the filter
coeficients, then you get the filtered signal just by multiplying boths.
