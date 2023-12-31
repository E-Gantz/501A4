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

//NOTE: according to http://soundfile.sapp.org/doc/WaveFormat/ " 16-bit samples are stored as 2's-complement signed integers, ranging from -32768 to 32767. "

//convert a 16-bit 2's-complement signed integer to a double scaled to -1.0 to +1.0
double sampleTodouble(int16_t value){
   int16_t maxVal = 32767;

   if (value < 0){
        double absValue = abs(value);
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
    if (value > 1.0){
        result = maxVal;
    }
    else if (value < -1.0){
        result = -maxVal;
    }
    else {
        result = rint((double)value * (double)maxVal);
    }

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

//taken from class notes
void convolve(std::vector<double> input, std::vector<double> filter, std::vector<double> &y){
    int N = input.size();
    int M = filter.size();
    int P = N+(M-1);
    /* Outer loop: process each input value x[n] in turn */
    for (int n = 0; n < N; n++) {
    /* Inner loop: process x[n] with each sample of h[n] */
        for (int m = 0; m < M; m++){
            y[n+m] = y[n+m] + (input[n] * filter[m]);
        }
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
    std::vector<double> convolved(doubleSamples.size() + (doubleIR.size() - 1), 0.0);

    //convolution stuff
    convolve(doubleSamples, doubleIR, convolved);

    //doubles back to ints
    std::vector<int16_t> outputSamples = doublesToSamples(convolved);

    //Write the converted samples to a new WAV file
    writeWavFile(outputFile, outputSamples, 1);

    std::cout << "Reverbification complete. Output file: " << outputFile << std::endl;

    return 0;
}