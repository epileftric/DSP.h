/*
 * Samples.h
 *
 *  Created on: Aug 7, 2016
 *      Author: epileftric 
 * 
 *   This namespace defines and implements the template  samples_t<n,type> for signal processing. 
 *  The general idea of it's implementation is to simplify the process of dealing with samples 
 *  and filtering.  
 *   
 *  Just by setting two samples_t you can use one for the samples you read and another for the filter 
 *  coeficients, then you get the filtered signal just by multiplying boths. 
 *  
 *  
 *  Example: 
 *  
 *     samples_t<3,float> samples,  
 *              filter; 
 *    float  output; 
 *  
 *     filter << 0.25 
 *     filter << 0.50 
 *     filter << 0.25 
 *     
 *     while(true){ 
 *        samples << readADC(); 
 *       output = samples * filter; 
 *        set_output( output ); 
 *       delay(); 
 *     } 
 *   
 */ 

#ifndef SRC_DSP_H_
#define SRC_DSP_H_

#include <string.h>
#include <math.h>

namespace DSP {

	template<unsigned int N=3,typename T=float>
	struct samples_t {
	 public:
	 
	 
		//~ Return the lastest sample
		inline T& last(void){
			return s[N-1];
		}
		//~ Return the oldest sample
		inline T& first(void){
			return s[0];
		}
		
		//~ Allows to access a sample with it's time index: -1 -2 -3 -4 ... -N
		inline T& operator [] (int i) {
			return s[N+i];
		}
		
		//~ Adds a new sample and rotates the others
		inline void operator << (T f){
			memmove(s,s+1,(N-1)*sizeof(T));
			s[N-1] = f;
		}
		//~ Removes one sample from the vector and rotates the others
		inline void operator >> (T &f){
			f = s[N-1];
			memmove(s+1,s,(N-1)*sizeof(T));
		}

		//~ Cross product between two samples vectors
		inline samples_t<N,T> operator ^ (samples_t<N,T> k){
			samples_t<N,T> ret;
			for( unsigned int i = 0; i < N ; i++)
				ret.s[i] = s[i] * k.S()[i];
			return ret;
		}
		//~ Cross product between the current and another sample vector
		inline void operator ^= (samples_t<N,T> k){
			for( unsigned int i = 0; i < N ; i++)
				s[i] *= k.S()[i];
		}

		//~ Inner product, returns type T
		inline T operator * (samples_t<N,T> k){
			T f = 0.0;
			for( unsigned int i = 0; i < N ; i++)
				f += s[i] * k.S()[i];
			return f;
		}

		template<typename R>
		inline T operator * (samples_t<N,R> k ){
			R f = 0.0;
			for( unsigned int i = 0; i < N ; i++)
				f += s[i] * k.S()[i];
			return (T)f;
		}
		
		//~ Adds a sample to the vector, increasing its size by 1
		//~ There's no += operator since this kind of operation changes the type
		inline samples_t<N+1,T> operator + (T f){
			samples_t<N+1,T> ret;
			memcpy( ret.S(), s, N*sizeof(T));
			ret.S()[N] = f;

			return ret;
		}
		
		//~ Adds a sample vector to the vector, increasing its size by M
		//~ There's no += operator since this kind of operation changes the type
		template<unsigned int M>
		inline samples_t<N+M,T> operator + (samples_t<M,T> &a){
			samples_t<N+M,T> ret;
			memcpy( ret.S() + 0, this->S(), N*sizeof(T));
			memcpy( ret.S() + N, a.S(), M*sizeof(T));
			return ret;
		}
		
		
		//~ Copies a vector into another
		inline samples_t<N,T>& operator = (samples_t<N,T> v){
			memcpy( s, v.s, N*sizeof(T));
			return *this;
		}

		//~ Clears the whole vector
		inline void clear(void){
			memset(s,0,N*sizeof(T));
		}
		
		//~ Returns the mean value of the samples in the vector
		inline T mean(void){
			return (T) sum()/N;
		}
		
		//~ Returns the total sum of the samples in the vector
		inline T sum(void){
			T a = (T) 0;
			for (unsigned int i = 0; i < N; i++)
				a += s[i];
			return a;
		}
		
		T min_sample(void){
			T m = s[0];
			for (size_t i = 0; i < N; i++) {
				m = s[i] < m ? s[i] : m;
			}
			return m;
		}

		T max_sample(void){
			T M = s[0];
			for (size_t i = 0; i < N; i++) {
				M = s[i] > M ? s[i] : M;
			}
			return M;
		}
		
		
		//~ Sorts the samples in an auxiliar vector and then returns the value that's in the midle
		T median(void){
			samples_t<N,T> sorted = this->sorted();
			return sorted.s[N/2];					
		}
		
		inline T stddev(void){
			T x_bar,
			  sd= (T) 0;

			x_bar = mean();
			for( int i = 0 ;  i < N ; i++ ){
				sd += pow(s[i] - x_bar, 2 );
			}
			sd /= (N-1);
			return sqrt(sd);
		}
		
		//~ Return the convolution from this vector with another of different size and any other type
		template<unsigned int M, typename T2>
		inline samples_t<N,T> conv(samples_t<M,T2> k){
			samples_t<N,T> ret;
			for (unsigned int i = 0; i < N; i++)	{
				T acum = (T)0;
				for (unsigned int j = 0; j < M; j++){
					acum += s[i + j] * k.S()[j];
				}
				ret <<  acum;
			}
			return ret;
		}
		
		
		//~ Sorts the current vector.
		inline void sort(void){
			quicksort(0,N-1);
		}
		
		//~ Returns a sorted copy of the current vector
		inline samples_t<N,T>& sorted(void){
			samples_t<N,T>& ret = *this;
			ret.quicksort(0,N-1);
			return ret;
		
		}
		
		
		//~ Used only for PID controllers, example:
		//! samples_t<3,float> errors,pid;
		//! float output = 0.0;
		//! pid.kPID(0.1,0.2,0.3);
		//! while(1){
			//! errors << ref - measurement;
			//! output +=  errors * pid;
			//! setOutput( output );
		//! }
		samples_t<3,T>& kPID(T kp, T ki, T kd){
			*this << (T)(kd); 			    //K2
			*this << (T)(-kp - (2 * kd) );  //K1
			*this << (T)(kp + ki + kd);     //K0

			return *this;
		}



		//! reference to the internal buffer
		inline T* S(void){
			return s;
		}
		
		private:
		
			T s[N]= {(T)0};

		void quicksort(int low, int high)	{
			int pivot, i, j;
			T temp;
			if(low < high) {
				pivot = low; // select a pivot element
				i = low;
				j = high;
				while(i < j) {
					// increment i till you get a number greater than the pivot element
					while(s[i] <= s[pivot] && i <= high)
					i++;
					// decrement j till you get a number less than the pivot element
					while(s[j] > s[pivot] && j >= low)
					j--;
					// if i < j swap the elements in locations i and j
					if(i < j) {
						temp = s[i];
						s[i] = s[j];
						s[j] = temp;
					}
				}

				// when i >= j it means the j-th position is the correct position
				// of the pivot element, hence swap the pivot element with the
				// element in the j-th position
				temp = s[j];
				s[j] = s[pivot];
				s[pivot] = temp;
				// Repeat quicksort for the two sub-arrays, one to the left of j
				// and one to the right of j
				quicksort(low, j-1);
				quicksort(j+1, high);
				}
			}
		};

	typedef samples_t<6,float> samples6f_t;
	typedef samples_t<5,float> samples5f_t;
	typedef samples_t<4,float> samples4f_t;
	typedef samples_t<3,float> samples3f_t;
	typedef samples_t<2,float> samples2f_t;

	typedef samples_t<3,float> pid_t;

};

#endif /* SRC_DSP_H_ */
