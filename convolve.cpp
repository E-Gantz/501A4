#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>
#include <cmath>
#include <string>

/*  Standard sample rate in Hz  */
#define SAMPLE_RATE       44100.0

/*  Standard sample size in bits  */
#define BITS_PER_SAMPLE   16

/*  Standard sample size in bytes  */		
#define BYTES_PER_SAMPLE  (BITS_PER_SAMPLE/8)

/*  Number of channels  */
#define MONOPHONIC        1
#define STEREOPHONIC      2

#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

//NOTE: according to http://soundfile.sapp.org/doc/WaveFormat/ " 16-bit samples are stored as 2's-complement signed integers, ranging from -32768 to 32767. "

//convert a 16-bit 2's-complement signed integer to a double scaled to -1.0 to +1.0
double sampleTodouble(int16_t value){
   int16_t maxVal = 32767;

   if (value < 0){
        int16_t absValue = abs(value);
        double result = (double)absValue / (double)maxVal;
        return -result;
   }
   else {
        double result = (double)value / (double)maxVal;
        return result;
   }
}

//convert a double scaled to -1.0 to +1.0 to a 16-bit 2's-complement signed integer
int16_t doubleToSample(double value){
    int16_t maxVal = 32767;
    int16_t result;
    // if (value > 1.0){
    //     result = maxVal;
    // }
    // else if (value < -1.0){
    //     result = -maxVal;
    // }
    // else {
    //     result = rint(value * maxVal);
    // }
    result = rint(value * maxVal);

    return result;
}

//convert the vector of wav samples to doubles
std::vector<double> samplesTodoubles(std::vector<int16_t> samples){
   std::vector<double> doubleSamples;
   for (int16_t sample : samples){
        double convertedSample = sampleTodouble(sample);
        doubleSamples.push_back(convertedSample);
   }
   return doubleSamples;
}

//convert the vector of doubles into wav samples
std::vector<int16_t> doublesToSamples(std::vector<double> doubleSamples){
   std::vector<int16_t> samples;
   for (double sample : doubleSamples){
        int16_t convertedSample = doubleToSample(sample);
        samples.push_back(convertedSample);
    }
   return samples;
}


std::vector<int16_t> readWavFile(char *filename){
    std::ifstream file(filename, std::ios::binary);

    if (!file.is_open()){
        fprintf(stdout, "File %s cannot be opened for reading\n", filename);
        return std::vector<int16_t>();
    }

    //read header
    char header[44];
    file.read(header, 44);

    //Check if the file is mono, 16-bit, 44.1 kHz
    //header[22] is numchannels, 34 is bits per sample, 24 and 26 is sample rate
    if (header[22] != 1 || header[34] != 16 || *reinterpret_cast<int*>(header + 24) != 44100){
        std::cerr << "Unsupported WAV file format. Must be mono, 16-bit, 44.1 kHz." << std::endl;
        return std::vector<int16_t>();
    }

    //store samples
    std::vector<int16_t> samples;
    int16_t sample;
    while (file.read(reinterpret_cast<char*>(&sample), sizeof(int16_t))){
        samples.push_back(sample);
    }

    file.close();
    return samples;
}

//Note: The following are from testtone.c given on D2L
size_t fwriteIntLSB(int data, FILE *stream)
{
    unsigned char array[4];

    array[3] = (unsigned char)((data >> 24) & 0xFF);
    array[2] = (unsigned char)((data >> 16) & 0xFF);
    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 4, stream);
}

size_t fwriteShortLSB(int16_t data, FILE *stream)
{
    unsigned char array[2];

    array[1] = (unsigned char)((data >> 8) & 0xFF);
    array[0] = (unsigned char)(data & 0xFF);
    return fwrite(array, sizeof(unsigned char), 2, stream);
}

void writeWaveFileHeader(int channels, int numberSamples,
                         double outputRate, FILE *outputFile)
{
    /*  Calculate the total number of bytes for the data chunk  */
    int dataChunkSize = channels * numberSamples * BYTES_PER_SAMPLE;
	
    /*  Calculate the total number of bytes for the form size  */
    int formSize = 36 + dataChunkSize;
	
    /*  Calculate the total number of bytes per frame  */
    int16_t frameSize = channels * BYTES_PER_SAMPLE;
	
    /*  Calculate the byte rate  */
    int bytesPerSecond = (int)ceil(outputRate * frameSize);

    /*  Write header to file  */
    /*  Form container identifier  */
    fputs("RIFF", outputFile);
      
    /*  Form size  */
    fwriteIntLSB(formSize, outputFile);
      
    /*  Form container type  */
    fputs("WAVE", outputFile);

    /*  Format chunk identifier (Note: space after 't' needed)  */
    fputs("fmt ", outputFile);
      
    /*  Format chunk size (fixed at 16 bytes)  */
    fwriteIntLSB(16, outputFile);

    /*  Compression code:  1 = PCM  */
    fwriteShortLSB(1, outputFile);

    /*  Number of channels  */
    fwriteShortLSB((int16_t)channels, outputFile);

    /*  Output Sample Rate  */
    fwriteIntLSB((int)outputRate, outputFile);

    /*  Bytes per second  */
    fwriteIntLSB(bytesPerSecond, outputFile);

    /*  Block alignment (frame size)  */
    fwriteShortLSB(frameSize, outputFile);

    /*  Bits per sample  */
    fwriteShortLSB(BITS_PER_SAMPLE, outputFile);

    /*  Sound Data chunk identifier  */
    fputs("data", outputFile);

    /*  Chunk size  */
    fwriteIntLSB(dataChunkSize, outputFile);
}

//Function to write audio data to a new WAV file
void writeWavFile(char *filename, std::vector<int16_t> samples, int numberOfChannels){
    FILE *outputFileStream = fopen(filename, "wb");
    if (outputFileStream == NULL){
        fprintf(stdout, "File %s cannot be opened for writing\n", filename);
        return;
    }

    int numberOfSamples = samples.size();

    /*  Write the WAVE file header  */
    writeWaveFileHeader(numberOfChannels, numberOfSamples,
                        SAMPLE_RATE, outputFileStream);

    //Write audio data (samples) to the file
    for (int i = 0; i < numberOfSamples; i++){
        fwriteShortLSB(samples[i], outputFileStream);
    }

    fclose(outputFileStream);
}

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
	if (j > i) {
	    SWAP(data[j], data[i]);
	    SWAP(data[j+1], data[i+1]);
	}
	m = nn;
	while (m >= 2 && j > m) {
	    j -= m;
	    m >>= 1;
	}
	j += m;
    }

    mmax = 2;
    while (n > mmax) {
	istep = mmax << 1;
	theta = isign * (6.28318530717959 / mmax);
	wtemp = sin(0.5 * theta);
	wpr = -2.0 * wtemp * wtemp;
	wpi = sin(theta);
	wr = 1.0;
	wi = 0.0;
	for (m = 1; m < mmax; m += 2) {
	    for (i = m; i <= n; i += istep) {
		j = i + mmax;
		tempr = wr * data[j] - wi * data[j+1];
		tempi = wr * data[j+1] + wi * data[j];
		data[j] = data[i] - tempr;
		data[j+1] = data[i+1] - tempi;
		data[i] += tempr;
		data[i+1] += tempi;
	    }
	    wr = (wtemp = wr) * wpr - wi * wpi + wr;
	    wi = wi * wpr + wtemp * wpi + wi;
	}
	mmax = istep;
    }
}

int nearestPower(int M){
    int i = 0;
    while (true){
        if (M <= pow(2, i)){
            return pow(2,i);
        }
        else {
            i++;
        }
    }
}

//taken from 311-318 in the smith text
void convolve(std::vector<double> input, std::vector<double> filter, std::vector<double> &y){
    int N = input.size();
    int M = filter.size();
    int outSize = nearestPower(M);
    int numChunkSamples = outSize-(M-1);
    int numChunks = N/numChunkSamples;//fix this
    int fftSize = outSize/2;

    double* XX = new double[outSize];
    std::vector<double> REX(fftSize, 0.0);
    std::vector<double> IMX(fftSize, 0.0);
    std::vector<double> REFR(fftSize, 0.0);
    std::vector<double> IMFR(fftSize, 0.0);
    std::vector<double> OLAP(M-1, 0.0);

    //use FFT to fill in REFR and IMFR
    // "After defining and initializing all
    // the arrays (lines 130 to 230), the first step is to calculate and store the
    // frequency response of the filter (lines 250 to 310). Line 260 calls a
    // mythical subroutine that loads the filter kernel into XX[0] through
    // XX[399], and sets XX[400] through XX[1023] to a value of zero."
    for (int i=0; i<M; i++){
        XX[i] = filter[i];
    }
    for (int i=M; i<outSize; i++){
        XX[i] = 0.0;
    }
    // The subroutine in line 270 is the FFT, transforming the 1024 samples held in
    // XX[ ] into the 513 samples held in REFR[ ] & IMFR[ ] 
    four1(XX-1, fftSize, 1);
    int counter = 0;
    for (int i=0; i<outSize-1; i+=2){
        REFR[counter] = XX[i];
        IMFR[counter] = XX[i+1];
        counter++;
    }

    for (int i=0; i<numChunks; i++){
        //use FFT to fill REX and IMX
        int inputOffset = i*numChunkSamples;
        for (int j=0; j<numChunkSamples; j++){
            XX[j] = input[inputOffset + j];
        }
        for (int j=numChunkSamples; j<outSize; j++){
            XX[j] = 0.0;
        }
        four1(XX-1, fftSize, 1);
        counter = 0;
        for (int j=0; j<outSize-1; j+=2){
            REX[counter] = XX[j];
            IMX[counter] = XX[j+1];
            counter++;
        }

        for (int j=0; j<fftSize; j++){
            double temp = (REX[j]*REFR[j]) - (IMX[j]*IMFR[j]);
            IMX[j] = (REX[j]*IMFR[j]) + (IMX[j]*REFR[j]);
            REX[j] = temp;
        }

        //use IFFT to combine REX and IMX into XX
        counter = 0;
        for (int j=0; j<fftSize; j++){
            XX[counter] = REX[i];
            counter ++;
            XX[counter] = IMX[i];
            counter++;
        }
        four1(XX-1, fftSize, -1);

        //Scale the output by N
        for (int j=0; j<outSize; j++){
            XX[j] = (double)XX[j] / (double)N;
        }

        for (int j=0; j<M-1; j++){
            XX[j] = XX[j] + OLAP[j];
        }

        for (int j = numChunkSamples; j<outSize; j++){
            OLAP[j-numChunkSamples] = XX[j];
        }

        for (int j=0; j<numChunkSamples; j++){
            y[inputOffset + j] = XX[j];
        }
    }

    for (int i=0; i<M-1; i++){
        y[N-1 + i] = OLAP[i];
    }
}

int main(int argc, char *argv[]){
    char *inputFile = NULL;
    char *IRFile = NULL;
    char *outputFile = NULL;

    if (argc != 4){
        std::cout << "Usage: " << argv[0] << " inputfile IRfile outputfile" << std::endl;
        return 1;
    }

    inputFile = argv[1];
    IRFile = argv[2];
    outputFile = argv[3];

    //get the data section of the input file as a vector
    std::vector<int16_t> wavSamples = readWavFile(inputFile);
    if (wavSamples.empty()){
        return 1;
    }
    std::vector<int16_t> irSamples = readWavFile(IRFile);
    if (irSamples.empty()){
        return 1;
    }

    //ints to doubles
    std::vector<double> doubleSamples = samplesTodoubles(wavSamples);
    std::vector<double> doubleIR = samplesTodoubles(irSamples);
    std::vector<double> convolved(doubleSamples.size() + doubleIR.size() - 1, 0.0);

    //convolution stuff
    convolve(doubleSamples, doubleIR, convolved);
    //doubles back to ints
    std::vector<int16_t> outputSamples = doublesToSamples(convolved);
    //std::vector<int16_t> outputSamples = doublesToSamples(doubleSamples);

    //Write the converted samples to a new WAV file
    writeWavFile(outputFile, outputSamples, 1);

    std::cout << "Reverbification complete. Output file: " << outputFile << std::endl;

    return 0;
}